"""What travels between blocks.

Six types carry everything the two arms do. They are immutable, so a corpus
that crosses five blocks leaves five comparable states rather than one that has
been edited in place five times, and an article can never lose a field on the
way to the model without that showing in a diff.

The dict shape the rest of the codebase uses is preserved exactly by from_dict
and as_dict, in both directions, because runs recorded on disk are read back
years later and a rename here would strand them.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable, Mapping, Optional

# Every field an article can carry. A block declares which of these it reads,
# and the composer checks they are there by then: `abstract` only exists once
# the hydrate step has run.
FIELDS = frozenset({"pmid", "title", "abstract", "found_by", "date", "tier"})
# What PubMed gives before abstracts are fetched.
SEARCHED = frozenset({"pmid", "title", "found_by", "date", "tier"})

VERDICTS = ("Novel", "Hypothesized", "Existing", "Established")


@dataclass(frozen=True)
class Query:
    """The question. Everything else in a run is derived from these two."""
    gene: str
    phenotype: str

    def __str__(self):
        return f"{self.gene} / {self.phenotype}"


@dataclass(frozen=True)
class Article:
    pmid: str
    title: str = ""
    abstract: str = ""
    found_by: tuple = ()
    date: str = ""
    tier: str = ""

    @classmethod
    def from_dict(cls, d):
        return cls(pmid=str(d["pmid"]), title=d.get("title") or "",
                   abstract=d.get("abstract") or "",
                   found_by=tuple(d.get("found_by") or ()),
                   date=d.get("pubdate") or "", tier=d.get("query_tier") or "")

    def as_dict(self):
        return {"pmid": self.pmid, "title": self.title,
                "abstract": self.abstract, "found_by": list(self.found_by),
                "pubdate": self.date, "query_tier": self.tier}


@dataclass(frozen=True)
class Corpus:
    """Articles in a declared order, unique by PMID.

    The order is part of the value, not an accident of which search answered
    first: it ranks the queries, so selection can break a tie on how an article
    was reached rather than on scheduling.

    `provided` is the set of fields the articles are known to carry. It is
    carried rather than guessed from the contents, since an empty abstract is a
    fact about one article and not about the corpus.
    """
    articles: tuple = ()
    provided: frozenset = SEARCHED

    @classmethod
    def of(cls, dicts, provided=SEARCHED):
        seen, out = set(), []
        for d in dicts:
            a = d if isinstance(d, Article) else Article.from_dict(d)
            if a.pmid not in seen:
                seen.add(a.pmid)
                out.append(a)
        return cls(tuple(out), frozenset(provided))

    def as_dicts(self):
        return [a.as_dict() for a in self.articles]

    @property
    def pmids(self):
        return tuple(a.pmid for a in self.articles)

    def __len__(self):
        return len(self.articles)

    def __iter__(self):
        return iter(self.articles)

    def __getitem__(self, i):
        return self.articles[i]

    def merged(self, other):
        """This corpus then the other, first sighting of a PMID winning.

        A field is only provided afterwards if both sides had it, which is what
        stops a hydrated corpus merged with a fresh search from claiming
        abstracts it half has.
        """
        return Corpus.of(list(self.articles) + list(other.articles),
                         self.provided & other.provided)

    def keep(self, pmids):
        keep = set(pmids)
        return replace(self, articles=tuple(a for a in self.articles
                                            if a.pmid in keep))

    def with_abstracts(self, by_pmid):
        return Corpus(tuple(replace(a, abstract=by_pmid.get(a.pmid, a.abstract))
                            for a in self.articles),
                      self.provided | {"abstract"})


@dataclass(frozen=True)
class Label:
    """What a grader said about one article."""
    relevance: Optional[int] = None
    kind: Optional[str] = None
    dropped: bool = False

    def as_tier(self):
        """The (relevance, kind) pair the reading order is written in."""
        return (self.relevance, self.kind)


@dataclass(frozen=True)
class Labels:
    """Grades keyed by PMID, kept beside the corpus and never inside it."""
    by_pmid: Mapping = field(default_factory=dict)

    @classmethod
    def from_tiers(cls, tiers):
        return cls({p: Label(relevance=t[0], kind=t[1])
                    for p, t in tiers.items()})

    @classmethod
    def from_dropped(cls, pmids):
        return cls({str(p): Label(dropped=True) for p in pmids})

    def tiers(self):
        return {p: lab.as_tier() for p, lab in self.by_pmid.items()
                if lab.relevance is not None}

    def dropped(self):
        return {p for p, lab in self.by_pmid.items() if lab.dropped}

    def merged(self, other):
        out = dict(self.by_pmid)
        for pmid, lab in other.by_pmid.items():
            was = out.get(pmid)
            out[pmid] = Label(
                relevance=lab.relevance if lab.relevance is not None
                else (was.relevance if was else None),
                kind=lab.kind or (was.kind if was else None),
                dropped=lab.dropped or bool(was and was.dropped))
        return Labels(out)

    def __len__(self):
        return len(self.by_pmid)


@dataclass(frozen=True)
class Distribution:
    """A hundred points spread over the four grades, as the reader returns them."""
    established: int = 0
    existing: int = 0
    hypothesized: int = 0
    novel: int = 0

    def as_dict(self):
        return {"Established": self.established, "Existing": self.existing,
                "Hypothesized": self.hypothesized, "Novel": self.novel}

    def top(self):
        d = self.as_dict()
        return max(d, key=d.get)


@dataclass(frozen=True)
class Verdict:
    call: str
    distribution: Optional[Distribution] = None
    justification: str = ""
    pmids: tuple = ()
    detail: Mapping = field(default_factory=dict)


@dataclass(frozen=True)
class Step:
    """What one block did, kept so a run can be read back without rerunning it."""
    block: str
    fingerprint: str
    model: Optional[str] = None
    seconds: float = 0.0
    note: str = ""
    input_tokens: int = 0
    output_tokens: int = 0
    usd: float = 0.0
    calls: int = 0


@dataclass(frozen=True)
class Trace:
    steps: tuple = ()

    def then(self, step):
        return Trace(self.steps + (step,))

    def as_list(self):
        return [vars(s) for s in self.steps]


def as_corpus(articles: Iterable, provided=SEARCHED):
    """Accept whatever the older code hands over, return a Corpus."""
    if isinstance(articles, Corpus):
        return articles
    return Corpus.of(articles, provided)
