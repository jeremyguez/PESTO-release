"""The one place a corpus becomes text for a model.

There were four ways of doing this: two flags on format_articles_for_llm, a
hand-built line for the title-only pass, and a corpus rebuilt with `"title": ""`
so the reader would not see titles. The last one was invisible in a diff and
changed what a published arm sent, which is the reason this file exists.

A block does not choose a format. It declares which fields it reads, and the
format follows from that.
"""
from __future__ import annotations

from ..services.pubmed_service import format_articles_for_llm
from .types import Corpus


def check_fields(corpus: Corpus, needs, block="block"):
    """Fail before the call rather than after it, with the reason.

    A block that asks for abstracts before the hydrate step gets a corpus of
    empty strings and a plausible, wrong answer. That is the failure worth
    making loud.
    """
    missing = frozenset(needs) - frozenset(corpus.provided)
    if missing:
        raise ValueError(
            f"{block} reads {sorted(missing)}, which the corpus does not carry "
            f"yet; it has {sorted(corpus.provided)}")


def render(corpus: Corpus, needs, check=True, block="block"):
    """Articles as the model will read them, and nothing more.

    Three shapes, picked by what is asked for rather than by a flag: title and
    provenance alone, the same with abstracts, or the bare `pmid title` lines
    the title-only pass uses.
    """
    needs = frozenset(needs)
    if check:
        check_fields(corpus, needs, block)
    if not len(corpus):
        return "No articles found." if "found_by" in needs else ""
    if "found_by" in needs:
        return format_articles_for_llm(corpus.as_dicts(),
                                       use_abstracts="abstract" in needs)
    return "\n".join(f"{a.pmid} {a.title}" for a in corpus)
