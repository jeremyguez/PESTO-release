"""The arms, written down.

An arm is a pipeline declared as data: the blocks it runs, in order, each with
its parameters and prompts. Its fingerprint hashes all of that, so a saved run
says exactly which pipeline produced it. The name is not hashed: renaming an
arm leaves every saved run valid.

Two families, each with a disease and a trait variant:

  abstracts        PESTO as in the paper, and the default. PubMed searched
                   under the phenotype's synonyms, every hit banded by Haiku
                   (strong human, weak human, non-human, no bearing), twenty
                   abstracts read by Opus, which spreads 100 points over Novel,
                   Hypothesized, Existing and Established.
  abstracts-trait  The same, with prompts that do not assume the phenotype is
                   a disease: Established also accepts replicated cohort and
                   GWAS evidence.
  fulltext         abstracts plus a Europe PMC search for articles that name
                   the gene only in their body, with the passages naming it
                   appended for the reader. Titles are sieved before abstracts
                   are fetched, so it costs about as much as abstracts.
  fulltext-trait   fulltext with the trait prompts of abstracts-trait.

By default the variant is chosen by one Opus call that asks whether the
phenotype is a disease. The other arms reproduce the paper's comparisons:

  titles           abstracts, titles only (PESTO-titles, Extended Data Fig. 2)
  knowledge        no search at all, the model's own knowledge (Fig. 2c,d)
  abstracts-bare   the reading prompt cut to four one-line definitions
  abstracts-open   four bins named by intensity only
  abstracts-score  one 0-100 novelty number instead of four bins
"""
from __future__ import annotations

from .spec import Arm
from .steps import (Anchors, Argmax, Bands, BandsGeneral, BareProbabilityRead,
                    BroadSearch, Broadened, ByStrength, ByStrengthThenFullText,
                    ContentWords, DedupeByTitle, DedupeByTitleNoNotices,
                    FetchAbstracts, FetchAbstractsAndExcerpts, FullTextSearch,
                    GeneNames, GeneticsNameSearch, KeepAll, KnowledgeRead,
                    NameSearch, NoveltyScoreRead, OpenProbabilityRead,
                    ProbabilityRead, ProbabilityReadGeneral, ScoreRule,
                    SkipAbstracts, Synonyms, SynonymsGeneral,
                    TitleKickSparingFullText, VariantSearch, WordSearch,
                    harness)

# What Opus is sent does not grow with the gene's fame. Haiku puts every
# article on one ladder (strong human, weak human, not human, not this gene),
# and a rule takes the whole strong band, drawn by the hash of the PMID if it
# overflows twenty, then fills from the weaker bands up to twenty. Band 1 is
# never read.
ABSTRACTS = Arm(
    name="abstracts",
    expand=(Synonyms("opus5"), Broadened("haiku"), ContentWords(),
            GeneNames(aliases=True)),
    sources=(NameSearch(primary=harness.NAME_HEAD_QUOTA), GeneticsNameSearch(),
             BroadSearch(), VariantSearch(), WordSearch(when="unnamed")),
    dedupe=DedupeByTitle(),
    hydrate=FetchAbstracts(),
    annotate=(Bands("haiku"),),
    select=ByStrength(budget=20, cap4=20),
    read=ProbabilityRead("opus5"),
    judge=Argmax(),
    note="PESTO in the paper: twenty abstracts read a pair",
)

# The disease wording taken out of the three calls that carry it, and nothing
# else moved: same searches, same bands, same budget.
ABSTRACTS_TRAIT = ABSTRACTS.but(
    name="abstracts-trait",
    expand=(SynonymsGeneral("opus5"), Broadened("haiku"), ContentWords(),
            GeneNames(aliases=True)),
    annotate=(BandsGeneral("haiku"),),
    read=ProbabilityReadGeneral("opus5"),
    note="abstracts, prompts that do not assume the phenotype is a disease",
)

# Abstracts never fetched: Haiku and Opus see titles only.
TITLES = ABSTRACTS.but(
    name="titles",
    hydrate=SkipAbstracts(),
    note="abstracts, titles throughout, abstracts never fetched",
)

# The three reading-prompt ablations of the Supplementary Note. Same articles
# as abstracts; only the instruction to the reader changes.
ABSTRACTS_BARE = ABSTRACTS.but(
    name="abstracts-bare",
    read=BareProbabilityRead("opus5"),
    note="abstracts, category definitions only",
)

ABSTRACTS_OPEN = ABSTRACTS.but(
    name="abstracts-open",
    read=OpenProbabilityRead("opus5"),
    note="abstracts, four bins by intensity only",
)

ABSTRACTS_SCORE = ABSTRACTS.but(
    name="abstracts-score",
    read=NoveltyScoreRead("opus5"),
    judge=ScoreRule(),
    note="abstracts, one 0-100 novelty number",
)

# Closed book: no search, no corpus.
KNOWLEDGE = Arm(
    name="knowledge",
    expand=(),
    sources=(),
    dedupe=DedupeByTitle(),
    hydrate=SkipAbstracts(),
    annotate=(),
    select=KeepAll(),
    read=KnowledgeRead("opus5"),
    judge=Argmax(),
    note="the model's own knowledge, no literature",
)

# abstracts plus one source and one thing fetched. PubMed cannot return an
# article whose abstract never names the gene, and the deciding paper is often
# one of those: the gene sits in a results table or one sentence of a
# discussion. FullTextSearch asks Europe PMC for the gene in the body and the
# phenotype in the title or abstract, and FetchAbstractsAndExcerpts appends the
# passages naming the gene so the article survives banding. A title sieve runs
# before abstracts are fetched, sparing the body-search hits, and the PubMed
# reading list goes from twenty to fourteen, so the six full-text places come
# at about the price of abstracts.
FULLTEXT = ABSTRACTS.but(
    name="fulltext",
    expand=ABSTRACTS.expand + (Anchors(),),
    sources=ABSTRACTS.sources + (FullTextSearch(),),
    dedupe=DedupeByTitleNoNotices(),
    hydrate=FetchAbstractsAndExcerpts(),
    annotate=(TitleKickSparingFullText("haiku"), Bands("haiku")),
    hydrate_at=1,
    select=ByStrengthThenFullText(budget=14, cap4=14, extra=6),
    note="abstracts plus the Europe PMC body search, fourteen plus six read",
)

FULLTEXT_TRAIT = FULLTEXT.but(
    name="fulltext-trait",
    expand=ABSTRACTS_TRAIT.expand + (Anchors(),),
    annotate=(TitleKickSparingFullText("haiku"), BandsGeneral("haiku")),
    read=ProbabilityReadGeneral("opus5"),
    note="fulltext, prompts that do not assume the phenotype is a disease",
)

ARMS = {arm.name: arm for arm in
        (ABSTRACTS, ABSTRACTS_TRAIT, TITLES, KNOWLEDGE,
         ABSTRACTS_BARE, ABSTRACTS_OPEN, ABSTRACTS_SCORE,
         FULLTEXT, FULLTEXT_TRAIT)}

# The names these arms carried while the paper was written. The result tables
# under paper/ record them in their `arm` column.
ALIASES = {
    "current": "abstracts",
    "current-general": "abstracts-trait",
    "prior": "knowledge",
    "current-bare": "abstracts-bare",
    "current-open": "abstracts-open",
    "current-score": "abstracts-score",
    "currentv2-cheap": "fulltext",
}

# A family is a disease arm and its trait variant. `auto` picks the variant.
FAMILIES = {"abstracts": ("abstracts", "abstracts-trait"),
            "fulltext": ("fulltext", "fulltext-trait")}
AUTO = "auto"


def get(name):
    name = ALIASES.get(name, name)
    if name not in ARMS:
        raise KeyError(f"no arm named {name!r}; there are {sorted(ARMS)}")
    return ARMS[name].validate()


def choose(phenotype, model="opus5", family="abstracts"):
    """The disease arm of a family if the name is a disease, its trait arm
    otherwise. One Opus call decides."""
    disease, trait = FAMILIES[family]
    kind = harness.phenotype_kind(phenotype, model)
    return ARMS[disease] if kind == "disease" else ARMS[trait]


def resolve(name, phenotype, model="opus5", family="abstracts"):
    """An explicit arm, or the variant of `family` the triage call picks."""
    if name in (None, "", AUTO):
        return choose(phenotype, model, family).validate()
    return get(name)
