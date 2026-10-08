"""The work the flow blocks call: query expansion, PubMed retrieval, the
annotation passes and the selection rules.

Every block in pesto.flow.steps wraps a function here, so a block's fingerprint
(its parameters and prompts) names what it does, while the body that does it
lives in this module. Changing a function body here changes what an arm reads
without moving its fingerprint: that is the one thing the fingerprint does not
protect.
"""
import hashlib
import os
import re

from .services.pubmed_service import search_pubmed_gene_phenotype  # noqa: F401
from .config import APP_ROOT
from .cost import thread_pool
from .services.llm_service import call_llm_with_usage
from .services.pubmed_service import format_articles_for_llm
from .utils.helpers import load_prompt

PROMPTS = os.path.join(APP_ROOT, "prompts")
BROAD_PROMPT = os.path.join(PROMPTS, "broad_term_prompt.txt")
# The aliases the arms are declared with, and the Anthropic models they mean.
# The aliases, not the model ids, enter the fingerprints.
MODELS = {"haiku": "claude-haiku-4-5", "opus": "claude-opus-5",
          "opus5": "claude-opus-5"}
SPECIFIC_QUOTA, BROAD_QUOTA, SECONDARY_QUOTA, OLDEST_QUOTA = 50, 30, 10, 10
DECOMPOSED_QUOTA = 10
# Unfiltered name search on a famous gene ranks cytokine biology above the
# association study that settled the question (IL6 / Kaposi, Foster 2000 at
# rank 181 of 262). AND-ing these words and keeping fifty puts that study
# in the corpus; the unfiltered head of twenty keeps the old case reports
# that never write them. Both numbers are parameters.
GENETICS_TERMS = ["polymorphism", "polymorphisms", "variant", "variants",
                  "mutation", "mutations"]
GENETICS_QUOTA = 50
NAME_HEAD_QUOTA = 20
# The disease name is the sharper signal but half the founding reports never
# write it, so pairing the gene with the words every variant report does use
# reaches them without opening the search to the gene's whole literature. On the
# 510 papers GenCC cites, this clause holds 71% of them in a fifth as many
# articles as the bare gene; the window is what the filter can afford to read.
VARIANT_TERMS = ["mutation", "mutations", "variant", "variants"]
VARIANT_QUOTA = 50
BATCH = 10
# Zero, because the rate limit is now held centrally in pubmed_service: one gate
# for the whole client, rather than each search sleeping on its own and getting
# it wrong in both directions at once.
DELAY_MS = 0
DISMISS = ("Dismiss papers that mention acronyms that are not the gene "
           "(homonym acronym).")


def slug(gene, phenotype):
    return f"{gene}_{re.sub(r'[^A-Za-z0-9]+', '_', phenotype)}".strip("_")


def ask(prompt, model, agent):
    return (call_llm_with_usage(prompt, MODELS.get(model, model), 0.0,
                                agent_name=agent) or {}).get("text", "")


def expand_synonyms(phenotype, model="opus5", prompt="phenotype_synonyms_prompt"):
    """The names a disease is written under, asked fresh every time.

    Both arms used to replay a list frozen elsewhere, the published one from the
    pair's first production run and the arithmetic one from a column of the
    benchmark table. That is right for a benchmark, whose search must not move
    between runs, and wrong for anything else: a pair that was never run in
    production and sits in no table searched its bare name, which is how DDX41
    and ERBB3 lost the papers that named them.

    One call rather than production's four. The four exist to pool 0.0, 0.1, 0.2
    and 0.3, and the Anthropic client drops temperature for every Opus from 4.x
    on, so on this model they are the same request four times over.
    """
    from .parsers.llm_parsers import parse_phenotype_synonyms_response
    text = ask(load_prompt(prompt).format(phenotype=phenotype),
               model, "phenotype_synonyms_agent")
    out, seen = [], set()
    for syn in parse_phenotype_synonyms_response(text):
        s = syn.strip()
        if s and s.lower() not in seen:
            seen.add(s.lower())
            out.append(s)
    return out


def broad_term(phenotype, model):
    tpl = open(BROAD_PROMPT, encoding="utf-8").read().strip()
    m = re.search(r"<term>(.*?)</term>",
                  ask(tpl.format(phenotype=phenotype), model, "broad_term"), re.S)
    term = (m.group(1).strip() if m else "")
    return "" if not term or term.upper() == "NONE" \
        or term.lower() == phenotype.lower() else term


# Words carrying no clinical content. Removing them turns a label no author
# writes, such as `complex neurodevelopmental disorder`, into the one word they
# do write. Numbers and inheritance modes go the same way.
FILLER = {"and", "or", "of", "the", "with", "without", "type", "autosomal",
          "recessive", "dominant", "linked", "complex", "disorder", "disease",
          "syndrome", "familial", "susceptibility", "related",
          # Modifiers that qualify a disease name without naming anything. As
          # alternatives in a query they match half of PubMed, so `primary`
          # would have carried the search on its own.
          "primary", "secondary", "congenital", "hereditary", "idiopathic",
          "severe", "mild", "early", "late", "onset", "infantile", "juvenile",
          "adult", "progressive", "chronic", "acute", "isolated"}


def decompose(phenotype, gene=None):
    """Content words of a label, for when the whole phrase matches nothing.

    Searched as alternatives rather than jointly: the two words of `visceral
    neuropathy` never appear in the same abstract of the papers that matter,
    since authors write about Hirschsprung disease instead.

    The gene symbol is stripped first. Labels of the form `B3GALT6-congenital
    disorder of glycosylation` otherwise split into `galt`, which is a different
    gene entirely.
    """
    text = phenotype.lower()
    if gene:
        text = re.sub(re.escape(gene.lower()), " ", text)
    toks = [t.strip("-") for t in re.split(r"[^A-Za-z-]+", text)]
    return [t for t in dict.fromkeys(toks) if len(t) > 3 and t not in FILLER]



# The bands of the cheap arm, one axis instead of two. What Opus is paid to do
# is weigh human evidence, so the only distinctions that change what it is sent
# are: is this human, is the disease the right one, is the article about this
# gene at all. Kind and relevance collapse into that single ladder.
BANDS = {4: "strong", 3: "weak", 2: "nonhuman", 1: "none"}
BAND_PROMPT = os.path.join(PROMPTS, "evidence_band_prompt.txt")
# How many articles Opus reads when the corpus has more than it needs. Every
# band 4 goes in whatever the count; the rest is filled to here and no further.
READING_BUDGET = 30


def phenotype_kind(phenotype, model="opus5"):
    """disease or other. Unparseable answers count as disease."""
    raw = ask(load_prompt("phenotype_kind_prompt").format(phenotype=phenotype),
              model, "phenotype_kind")
    words = re.findall(r"[A-Za-z]+", (raw or "").lower())
    if "other" in words and "disease" not in words:
        return "other"
    if words and words[0] == "other":
        return "other"
    if words and words[-1] == "other":
        return "other"
    return "disease"


def band(gene, phenotype, arts, model, batch=BATCH, synonyms=None, prompt=None):
    """Each article on one ladder: strong human, weak human, not human, not this.

    A grade here is not a description of the article, it is a place in the
    queue for Opus's attention, and the prompt says so: when two bands are
    open, take the higher one. Batches of fixed size, so the answer does not
    depend on how long the corpus is.
    """
    path = (os.path.join(PROMPTS, f"{prompt}.txt") if prompt else BAND_PROMPT)
    tpl = open(path, encoding="utf-8").read().strip()
    names = "\n".join(f"  {s}" for s in (synonyms or []))
    chunks = [arts[i:i + batch] for i in range(0, len(arts), batch)]

    def label(chunk):
        text = ask(tpl.format(
            gene_name=gene, phenotype=phenotype,
            synonym_note=names or f"  {phenotype}",
            articles_list=format_articles_for_llm(chunk, use_abstracts=True)),
            model, "harness_band")
        seen = {a["pmid"] for a in chunk}
        out = {}
        for block in re.findall(r"<article>(.*?)</article>", text, re.S):
            pmid = re.search(r"<pmid>\s*(\d+)\s*</pmid>", block)
            grade = re.search(r"<band>\s*([1-4])\s*</band>", block)
            if pmid and grade and pmid.group(1) in seen:
                out[pmid.group(1)] = int(grade.group(1))
        return out

    bands = {}
    if chunks:
        with thread_pool(max_workers=min(len(chunks), 6)) as pool:
            for part in pool.map(label, chunks):
                bands.update(part)
    # An article the model skipped is not evidence that it is worthless, and the
    # prompt's own rule applies to us as well: unsure means the higher band.
    for a in arts:
        bands.setdefault(a["pmid"], 2)
    return bands


def select_bands(arts, bands, budget=READING_BUDGET, cap4=None):
    """Every strong paper, then as much context as the budget still allows.

    Two knobs, both meant to be swept. `budget` is where the filling stops, and
    `cap4` is whether the strong band is itself capped: None reads all of it,
    however many there are, on the argument that a gene with sixty human reports
    is exactly the gene whose sixty reports settle the question.

    Band 1 is never read at any budget. It is the one claim the classifier was
    told to make only when certain, and honouring it is what lets the rest be
    generous.

    When the strong band overflows its cap, the survivors are drawn by the hash
    of their PMID. Taking the first twenty instead would take them in the order
    the searches happened to return, which is neither the order they were
    published in nor the order they matter in, and would quietly hand the
    reading list to whichever query ran first. The hash is arbitrary in the same
    way a draw is, and gives the same twenty every time the pair is run.
    """
    rank = {a["pmid"]: i for i, a in enumerate(arts)}
    by = {b: [a for a in arts if bands.get(a["pmid"], 2) == b] for b in (4, 3, 2)}
    chosen = list(by[4])
    if cap4 and len(chosen) > cap4:
        chosen = sorted(chosen,
                        key=lambda a: hashlib.sha256(
                            str(a["pmid"]).encode()).hexdigest())[:cap4]
    for lower in (3, 2):
        if len(chosen) >= budget:
            break
        chosen += by[lower][:budget - len(chosen)]
    return sorted(chosen, key=lambda a: rank[a["pmid"]])


KICK_PROMPT = os.path.join(PROMPTS, "title_kick_prompt.txt")


def kick_titles(gene, phenotype, arts, model, batch=BATCH):
    """PMIDs to remove on the strength of the title alone.

    A first stage before the grader, which is the expensive one because it reads
    abstracts. Grading from titles was tried instead and graded badly: the title
    says which gene and which disease an article is about but rarely whether
    anyone carried a variant, so two hundred articles the grader called direct
    human evidence slid one band down and the reading lists diverged.

    Eliminating is the half of the job a title can do. A title that names a
    different gene, or a subject unrelated to both halves of the question,
    settles the matter on its own, and 39% of the corpus goes that way.

    Asymmetric on purpose: what this stage removes, nothing downstream can bring
    back, so it is asked for certainty rather than for judgement.
    """
    tpl = open(KICK_PROMPT, encoding="utf-8").read().strip()
    chunks = [arts[i:i + batch] for i in range(0, len(arts), batch)]

    def judge(chunk):
        text = ask(tpl.format(
            gene_name=gene, phenotype=phenotype,
            articles_list="\n".join(f"{a['pmid']} {a.get('title', '')}"
                                    for a in chunk)), model, "harness_kick")
        block = re.search(r"<remove>(.*?)</remove>", text, re.S)
        seen = {a["pmid"] for a in chunk}
        # A batch names only its own articles; anything else is the model
        # recalling a PMID from elsewhere.
        return set(re.findall(r"\d{6,9}", block.group(1) if block else "")) & seen

    with thread_pool(max_workers=min(len(chunks), 6)) as pool:
        return set().union(*pool.map(judge, chunks)) if chunks else set()
