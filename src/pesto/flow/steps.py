"""The blocks themselves, each wrapping work that lives in pesto.harness.

A block adds a declaration of what it reads, of what prompts it sends, and of
the parameters that make it what it is. Those three make its fingerprint.
"""
from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

from .. import harness
from ..services import fulltext_service
from ..services.llm_service import call_llm_with_usage
from ..services.pubmed_service import (deduplicate_articles, fetch_abstracts,
                                       publication_notices)
from ..utils.data_loader import data_loader
from ..utils.helpers import load_prompt
from .render import render
from .spec import Block, Judge, Reader, optional
from .types import Corpus, Distribution, Labels, Verdict

# ---------------------------------------------------------------- expansion
# Expanders fill in the names a search will be run under. They write into a
# plain dict so a source can ask for the one it needs by name.


@dataclass(frozen=True)
class Synonyms(Block):
    """Names the disease is written under.

    `model=None` replays a list recorded elsewhere, which is what a benchmark
    needs: its search must not move between runs. Anything else expands live,
    because a pair that was never run in production and sits in no table would
    otherwise be searched on its bare name, and that is how DDX41 and ERBB3 lost
    the papers that named them.
    """
    model: str = None
    uses = ("phenotype_synonyms_prompt",)

    def run(self, query, terms, given=None):
        if self.model is None:
            return {"synonyms": tuple(given or ())}
        return {"synonyms": tuple(harness.expand_synonyms(
            query.phenotype, self.model, prompt=self.uses[0]))}

    def fingerprint(self):
        # A replayed list is not this prompt's output, so the prompt is not part
        # of the block's identity when it is not sent.
        from .spec import digest, prompt_sha
        prompts = {p: prompt_sha(p) for p in self.uses} if self.model else {}
        return digest({"kind": self.kind, "params": self.params(),
                       "prompts": prompts})


@dataclass(frozen=True)
class Broadened(Block):
    """The disease name with its subtype index stripped, or nothing."""
    model: str = "haiku"
    uses = ("broad_term_prompt",)

    def run(self, query, terms, given=None):
        return {"broad": harness.broad_term(query.phenotype, self.model)}


@dataclass(frozen=True)
class Anchors(Block):
    """One or two words authors write for the phenotype, for the full-text pass.

    Synonyms are exact phrases, and a phrase must appear word for word. An
    article that describes the phenotype in its own words matches none of the
    generated phrases and never reaches the corpus, even when Europe PMC
    indexes the gene in its body; one or two plain words find it.
    """
    model: str = "haiku"
    most: int = 4
    uses = ("phenotype_anchors_prompt",)

    def run(self, query, terms, given=None):
        text = harness.ask(load_prompt(self.uses[0]).format(
            phenotype=query.phenotype), self.model, "phenotype_anchors")
        # An anchor naming the gene turns "gene in the body, phenotype in the
        # abstract" into "gene anywhere": MC1R came back as an anchor for red
        # hair and took the query from 11 hits to 928.
        genes = {n.lower() for n in (query.gene,) + tuple(terms.get("gene_names") or ())}
        out, seen = [], set()
        for a in re.findall(r"<anchor>(.*?)</anchor>", text or "", re.S):
            a = a.strip()
            if any(w.lower() in genes for w in re.findall(r"[\w-]+", a)):
                continue
            if a and len(a.split()) <= 2 and a.lower() not in seen:
                seen.add(a.lower())
                out.append(a)
        return {"anchors": tuple(out[:self.most])}


@dataclass(frozen=True)
class ContentWords(Block):
    """The label's words, minus the gene symbol and the clinically empty ones."""

    def run(self, query, terms, given=None):
        return {"words": tuple(harness.decompose(query.phenotype, query.gene))}


@dataclass(frozen=True)
class GeneNames(Block):
    """Symbols and approved names the gene is published under.

    Authors write dystrophin, not DMD, so this is not a small difference, which
    is exactly why it is declared rather than implied.
    """
    aliases: bool = True

    def run(self, query, terms, given=None):
        return {"gene_names": tuple(data_loader.get_search_names_for_gene(query.gene))
                if self.aliases else ()}


# ------------------------------------------------------------------ sources
# One block, one PubMed pass, one quota. The list of them is the search, and
# the order of the list is the order hits are merged in, so an article's rank
# records which query reached it first.


@dataclass(frozen=True)
class Search(Block):
    """A pass over PubMed. Subclasses only say which terms they run under."""
    primary: int = 0
    secondary: int = 0
    oldest: int = 0
    tag: str = ""

    # A pass that only looks at the question can be run before any answer comes
    # back; one that looks at what the other passes found has to wait for them.
    deferred = False

    def terms_for(self, query, terms):
        raise NotImplementedError

    def applies(self, corpus, terms):
        return True

    def run(self, query, terms):
        words = self.terms_for(query, terms)
        if not words:
            return Corpus()
        found, _ = harness.search_pubmed_gene_phenotype(
            query.gene, query.phenotype, list(words),
            max_results_primary=self.primary,
            max_results_secondary=self.secondary,
            max_results_oldest=self.oldest,
            delay_ms=harness.DELAY_MS,
            gene_aliases=list(terms.get("gene_names") or ()))
        if self.tag:
            for a in found:
                a["query_tier"] = self.tag
        return Corpus.of(found)


@dataclass(frozen=True)
class NameSearch(Search):
    """The disease under every name it is written by. The sharpest query there
    is, and the one half the founding reports never match."""
    primary: int = harness.SPECIFIC_QUOTA
    secondary: int = harness.SECONDARY_QUOTA
    oldest: int = harness.OLDEST_QUOTA

    def terms_for(self, query, terms):
        return (query.phenotype,) + tuple(terms.get("synonyms") or ())


@dataclass(frozen=True)
class GeneticsNameSearch(Search):
    """The name search restricted to papers that use the vocabulary of
    genetic evidence.

    The unfiltered name pass ranks cytokine biology above a 2000 association
    study (IL6 / Kaposi, Foster 2000, rank 181 of 262). AND-ing these words
    and keeping fifty puts that study in the corpus. The unfiltered head is
    a separate NameSearch, so old case reports that never write mutation
    still arrive."""
    primary: int = harness.GENETICS_QUOTA
    tag: str = "genetics"

    def terms_for(self, query, terms):
        return (query.phenotype,) + tuple(terms.get("synonyms") or ())

    def run(self, query, terms):
        words = self.terms_for(query, terms)
        if not words:
            return Corpus()
        found, _ = harness.search_pubmed_gene_phenotype(
            query.gene, query.phenotype, list(words),
            max_results_primary=self.primary,
            max_results_secondary=self.secondary,
            max_results_oldest=self.oldest,
            delay_ms=harness.DELAY_MS,
            gene_aliases=list(terms.get("gene_names") or ()),
            require_terms=list(harness.GENETICS_TERMS))
        if self.tag:
            for a in found:
                a["query_tier"] = self.tag
        return Corpus.of(found)


@dataclass(frozen=True)
class BroadSearch(Search):
    """The same disease with its subtype index dropped."""
    primary: int = harness.BROAD_QUOTA
    oldest: int = harness.OLDEST_QUOTA

    def terms_for(self, query, terms):
        broad = terms.get("broad")
        return (broad,) if broad else ()


@dataclass(frozen=True)
class FullTextSearch(Block):
    """The gene in the body, the disease in the title or abstract.

    The one pass that does not go through PubMed, and the only one that can
    reach an article whose abstract never names the gene. That class is not
    marginal: a genome-wide association study often reports a gene only in a
    results table or a supplement. Every PubMed pass here, however it is
    worded, is blind to them.

    It runs on the synonyms the expander already produced and on the anchor
    words, so the disease is matched under the words authors actually use and
    not only under the name the question was asked in.
    """

    limit: int = 25
    aliases: int = 4
    deferred = False

    def applies(self, corpus, terms):
        return True

    def run(self, query, terms):
        found = fulltext_service.search_fulltext(
            query.gene, list(terms.get("gene_names") or ())[:self.aliases],
            query.phenotype, list(terms.get("synonyms") or ()),
            limit=self.limit, anchors=list(terms.get("anchors") or ()))
        return Corpus.of(found)


@dataclass(frozen=True)
class VariantSearch(Search):
    """The gene beside the vocabulary of variant reports.

    The one pass that does not depend on the disease being named correctly. On
    the 510 papers GenCC cites it holds 71% of them, where the disease name
    alone reached 42%.
    """
    primary: int = harness.VARIANT_QUOTA
    tag: str = "variant"

    def terms_for(self, query, terms):
        return tuple(harness.VARIANT_TERMS)


@dataclass(frozen=True)
class WordSearch(Search):
    """The label's content words as alternatives, for names no author writes.

    `when="unnamed"` fires only if nothing matched the disease name itself,
    which is the case it was written for. Run for every pair instead, it also
    fires where the exact search worked, and the loose articles it adds push the
    middle of the GenCC scale towards Established.
    """
    primary: int = harness.DECOMPOSED_QUOTA
    oldest: int = harness.DECOMPOSED_QUOTA
    when: str = "unnamed"

    @property
    def deferred(self):
        return self.when != "always"

    def terms_for(self, query, terms):
        return tuple(terms.get("words") or ())

    def applies(self, corpus, terms):
        if self.when == "always":
            return True
        # `primary` means an article matched the disease name rather than the
        # gene beside a generic word. The variant pass never looks at the name,
        # so its hits say nothing about whether the name was matched.
        return not any(a.tier == "primary" for a in corpus)


# ------------------------------------------------------------- corpus shape


@dataclass(frozen=True)
class DedupeByTitle(Block):
    """Two PMIDs for one paper, the later one kept. A preprint and its journal
    version are one observation, and counting them twice is replication."""

    def run(self, corpus):
        return Corpus.of(deduplicate_articles(corpus.as_dicts()), corpus.provided)


@dataclass(frozen=True)
class DedupeByTitleNoNotices(DedupeByTitle):
    """The same, except that where two records share a title an erratum or a
    retraction notice never wins. Keeping the later record kept the AJHG
    erratum of the EIF4A2 paper (36868207) over the paper itself (36528028),
    so the reader got the title and an empty abstract. Only colliding titles
    are looked up, so a corpus without a collision costs no request."""

    def run(self, corpus):
        articles = corpus.as_dicts()
        by_title = {}
        for a in articles:
            by_title.setdefault((a.get("title") or "").strip().lower(), []).append(a)
        colliding = [str(a["pmid"]) for group in by_title.values() if len(group) > 1
                     for a in group]
        notices = publication_notices(colliding) if colliding else set()
        kept = [a for a in articles if str(a["pmid"]) not in notices]
        return Corpus.of(deduplicate_articles(kept), corpus.provided)


@dataclass(frozen=True)
class FetchAbstracts(Block):
    """Abstracts from PubMed, which is the only step that adds a field.

    Where it sits is a real choice. The title pass reads no abstracts, so
    placing this after it fetches only for the articles that survived.

    `sections` is in the fingerprint on purpose. The fingerprint hashes
    parameters and prompts, not the body of fetch_abstracts, so a change in
    how much of a structured abstract that function returns would otherwise
    leave the run cache serving corpora read the other way. "all" joins every
    section.
    """

    sections: str = "all"

    def run(self, corpus):
        return corpus.with_abstracts(fetch_abstracts(list(corpus.pmids)))


@dataclass(frozen=True)
class FetchAbstractsAndExcerpts(FetchAbstracts):
    """Abstracts, and for the body-only articles the sentences naming the gene.

    A full-text hit usually has an abstract that never says the gene, which is
    why PubMed could not find it. Sent on as it stands, Bands reads an abstract
    about a cohort with no gene in it and puts the article in band 1, and
    the arm pays for a search whose results it then throws away. So the passages
    around the gene are appended to the abstract, marked as full text, and the
    reader is told where they came from.

    Only articles the full-text pass found are fetched, and only when the gene
    is not already in title or abstract: the excerpt is there to repair what
    PubMed could not show, not to enlarge what it could.
    """

    window: int = 420
    most: int = 3
    wants_query = True

    def run(self, corpus, query=None):
        corpus = super().run(corpus)
        if query is None:
            return corpus
        names = [query.gene] + list(
            data_loader.get_search_names_for_gene(query.gene) or ())
        needle = re.compile(rf"\b{re.escape(query.gene)}\b", re.I)
        wanted = [a.pmid for a in corpus
                  if a.tier == "fulltext"
                  and not needle.search(f"{a.title} {a.abstract}")]
        if not wanted:
            return corpus
        pmcids = fulltext_service.pmcids_for(wanted)
        # One article at a time this was the slowest step of the arm, a median
        # of 31 s a pair and up to 141 s: up to 25 full texts, each allowed 45 s.
        with ThreadPoolExecutor(max_workers=min(len(wanted), 8)) as pool:
            found = pool.map(lambda pmid: fulltext_service.gene_excerpts(
                pmcids.get(pmid, ""), query.gene, names,
                window=self.window, most=self.most), wanted)
            added = {pmid: " ".join(pieces)
                     for pmid, pieces in zip(wanted, found) if pieces}
        if not added:
            return corpus
        merged = {}
        for a in corpus:
            if a.pmid in added:
                head = (a.abstract + " ") if a.abstract else ""
                merged[a.pmid] = (f"{head}[full text, passages naming "
                                  f"{query.gene}] {added[a.pmid]}")
        return corpus.with_abstracts({**{a.pmid: a.abstract for a in corpus},
                                      **merged})


@dataclass(frozen=True)
class SkipAbstracts(Block):
    """Mark abstracts present and empty. Nothing is fetched, and the formatter
    omits an empty abstract, so every later block sees titles only.
    """

    def run(self, corpus):
        return corpus.with_abstracts({})


# --------------------------------------------------------------- annotators
# An annotator returns grades keyed by PMID and never touches the corpus. What
# it marks dropped is removed by the runner before the next annotator, so a
# title pass narrows what the abstract pass has to pay for.


@dataclass(frozen=True)
class TitleKick(Block):
    """Removal on the title alone, asked for certainty rather than judgement."""
    model: str = "haiku"
    batch: int = harness.BATCH
    needs = frozenset({"pmid", "title"})
    uses = ("title_kick_prompt",)

    def run(self, query, corpus, terms=None):
        return Labels.from_dropped(harness.kick_titles(
            query.gene, query.phenotype, corpus.as_dicts(), self.model,
            self.batch))


@dataclass(frozen=True)
class TitleKickSparingFullText(TitleKick):
    """The title sieve, never shown the articles the body search found.

    The sieve removes a title that names a different gene, and a body-only hit
    is by definition an article whose title and abstract do not name this
    one, such as a genome-wide association study that lists it in a table.
    Shown those, it would remove exactly what the full-text pass was added for.
    """

    def run(self, query, corpus, terms=None):
        body = {a.pmid for a in corpus if a.tier == "fulltext"}
        rest = corpus.keep(p for p in corpus.pmids if p not in body)
        if not len(rest):
            return Labels()
        return super().run(query, rest, terms)


@dataclass(frozen=True)
class Bands(Block):
    """One ladder instead of two: strong human, weak human, not human, not this.

    Written for the arm that pays Opus by the article. The grader describes an
    article; this one places it in a queue, and says so in its prompt, so that
    an article it cannot make out is pushed up rather than dropped.
    """
    model: str = "haiku"
    batch: int = harness.BATCH
    needs = frozenset({"pmid", "title", "found_by", "abstract"})
    uses = ("evidence_band_prompt",)

    def run(self, query, corpus, terms=None):
        bands = harness.band(query.gene, query.phenotype, corpus.as_dicts(),
                             self.model, self.batch,
                             synonyms=(terms or {}).get("synonyms"),
                             prompt=self.uses[0])
        return Labels.from_tiers({p: (b, harness.BANDS[b])
                                  for p, b in bands.items()})


# ---------------------------------------------------------------- selectors


@dataclass(frozen=True)
class KeepAll(Block):
    """Whatever survived the annotators is what gets read."""

    def run(self, corpus, labels):
        return corpus


@dataclass(frozen=True)
class ByStrength(Block):
    """All the strong human evidence, then context up to a budget and no further.

    The cap it puts on a large corpus is the point: what Opus is sent stops
    depending on how much has been written about the gene. `cap4=None` leaves
    the strong band uncapped, since sixty human reports are what settles a
    question rather than what pads a bill.
    """
    budget: int = harness.READING_BUDGET
    cap4: int = None

    def run(self, corpus, labels):
        bands = {p: t[0] for p, t in labels.tiers().items()}
        kept = harness.select_bands(corpus.as_dicts(), bands, self.budget,
                                    self.cap4)
        return corpus.keep(a["pmid"] for a in kept)


@dataclass(frozen=True)
class ByStrengthThenFullText(ByStrength):
    """The PubMed selection first, then the body-only articles beside it.

    Adding a source without adding room is how a better search makes a worse
    answer. A first full-text run measured it: when the strong band
    overflowed, the hash draw inside `select_bands` evicted the article that
    carried the verdict, and a pair `abstracts` called Existing came back
    Novel. The full-text hits were not wrong; they simply took the seats.

    So the budget is spent on the PubMed corpus first, exactly as a ByStrength
    selection would, and the body-only articles are appended afterwards up to
    `extra`, strongest band first.
    """

    extra: int = 6

    def run(self, corpus, labels):
        body = {a.pmid for a in corpus if a.tier == "fulltext"}
        base = super().run(Corpus(tuple(a for a in corpus if a.pmid not in body),
                                  corpus.provided), labels)
        bands = {p: t[0] for p, t in labels.tiers().items()}
        add = sorted((a for a in corpus if a.pmid in body),
                     key=lambda a: -bands.get(a.pmid, 0))[:self.extra]
        return corpus.keep([a.pmid for a in base] + [a.pmid for a in add])


# ------------------------------------------------------------------ readers


@dataclass(frozen=True)
class ProbabilityRead(Reader):
    """One call over the whole corpus, returning a hundred points spread over
    the four grades.

    The prompt is the probabilities template with the closing rule of the
    assessment prompt appended, and the assertion below is why that is not
    left to memory.
    """
    model: str = "opus5"
    # Sentences put before the articles, one per line: --note, and the
    # instruction --max-date adds.
    notes: tuple = optional()
    evidence = "distribution"
    needs = frozenset({"pmid", "title", "found_by", "abstract"})
    uses = ("novelty_probabilities", "novelty_assessment_prompt")

    def template(self):
        head = load_prompt("novelty_probabilities").rstrip()
        _, _, rule = load_prompt("novelty_assessment_prompt").rpartition("\nIMPORTANT:")
        text = head + "\n\nIMPORTANT:" + rule.rstrip() + "\n"
        assert "at least two such PMIDs" in text, "closing rule not carried over"
        return text

    def run(self, query, corpus):
        resp = call_llm_with_usage(
            self.template().format(
                gene_name=query.gene, phenotype=query.phenotype,
                articles_list=render(corpus, self.needs, block=self.kind),
                alias_note="".join(n.strip() + "\n" for n in self.notes),
                dismiss_note=harness.DISMISS),
            harness.MODELS.get(self.model, self.model), 0.0,
            agent_name="gencc_probs") or {}
        return resp.get("text", "") or "", resp


@dataclass(frozen=True)
class BareProbabilityRead(ProbabilityRead):
    """Same reader as abstracts, with the category definitions only."""
    uses = ("novelty_probabilities_bare",)

    def template(self):
        return load_prompt("novelty_probabilities_bare")


@dataclass(frozen=True)
class OpenProbabilityRead(ProbabilityRead):
    """Four bins with intensity wording only: no species, no study count."""
    uses = ("novelty_probabilities_open",)

    def template(self):
        return load_prompt("novelty_probabilities_open")


@dataclass(frozen=True)
class NoveltyScoreRead(ProbabilityRead):
    """One novelty number from 0 to 100, no four-bin split."""
    evidence = "score"
    uses = ("novelty_score",)

    def template(self):
        return load_prompt("novelty_score")


@dataclass(frozen=True)
class SynonymsGeneral(Synonyms):
    """Same expansion as abstracts, without treating the name as a disease."""
    uses = ("phenotype_synonyms_general",)


@dataclass(frozen=True)
class BandsGeneral(Bands):
    """Same four bands as abstracts; the prompt no longer says 'disease'."""
    uses = ("evidence_band_general",)


@dataclass(frozen=True)
class ProbabilityReadGeneral(ProbabilityRead):
    """Same reader as abstracts; Established also accepts human cohorts."""
    uses = ("novelty_probabilities_general", "novelty_assessment_general")

    def template(self):
        head = load_prompt("novelty_probabilities_general").rstrip()
        _, _, rule = load_prompt("novelty_assessment_general").rpartition(
            "\nIMPORTANT:")
        text = head + "\n\nIMPORTANT:" + rule.rstrip() + "\n"
        assert "at least two such PMIDs" in text, "closing rule not carried over"
        return text


@dataclass(frozen=True)
class KnowledgeRead(Reader):
    """Closed book: the four grades from the model's knowledge, no corpus."""
    model: str = "opus5"
    evidence = "distribution"
    needs = frozenset()
    uses = ("knowledge_probabilities",)

    def template(self):
        return load_prompt("knowledge_probabilities")

    def run(self, query, corpus):
        resp = call_llm_with_usage(
            self.template().format(
                gene_name=query.gene, phenotype=query.phenotype),
            harness.MODELS.get(self.model, self.model), 0.0,
            agent_name="knowledge_prior") or {}
        return resp.get("text", "") or "", resp


# ------------------------------------------------------------------- judges
# No model, no network. A judgement can be replayed from a saved run and argued
# with.


@dataclass(frozen=True)
class Argmax(Judge):
    """The grade holding the most points."""
    reads = "distribution"

    def run(self, query, text, corpus=None):
        probs = {}
        for cat in ("Established", "Existing", "Hypothesized", "Novel"):
            m = re.search(rf"<{cat.lower()}>\s*(\d+)\s*</{cat.lower()}>", text)
            if m:
                probs[cat] = int(m.group(1))
        if not probs:
            return Verdict(call="PARSE_FAIL")
        dist = Distribution(established=probs.get("Established", 0),
                            existing=probs.get("Existing", 0),
                            hypothesized=probs.get("Hypothesized", 0),
                            novel=probs.get("Novel", 0))
        just = re.search(r"<justification>(.*?)</justification>", text, re.S)
        return Verdict(call=max(probs, key=probs.get), distribution=dist,
                       justification=(just.group(1).strip().replace("\t", " ")
                                      if just else ""),
                       detail=dist.as_dict())


@dataclass(frozen=True)
class ScoreRule(Judge):
    """The 0-100 novelty number put on the four labels.

    The same mapping as the mean of a distribution: 1 + 3 (100 - novelty) / 100,
    rounded and clipped to 1-4, so 0 is Established and 100 is Novel.
    """
    reads = "score"

    def run(self, query, text, corpus=None):
        m = re.search(r"<novelty>\s*([0-9]+(?:\.[0-9]+)?)\s*</novelty>",
                      text or "", re.I)
        if not m or not 0 <= float(m.group(1)) <= 100:
            return Verdict(call="PARSE_FAIL")
        novelty = float(m.group(1))
        level = min(4, max(1, round(1 + 3 * (100 - novelty) / 100)))
        return Verdict(call=("Novel", "Hypothesized", "Existing",
                             "Established")[level - 1],
                       detail={"novelty": novelty})
