"""Running an arm: the only place the order of the stages is written down.

The stages are fixed and the blocks are not. That is the trade this whole module
makes: what varies between arms varies inside a declaration that can be read at a
glance, and the sequence itself exists once, here, rather than in six scripts
that had each drifted a little from the others.
"""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

from .types import Corpus, Labels, Query, Step, Trace, Verdict

# An answer the model could not format is not one to trust for its numbers:
# the reading is asked again, and only the last answer is kept.
READ_TRIES = 3


@dataclass
class Result:
    query: Query
    arm: str
    fingerprint: str
    verdict: Verdict = None
    corpus: Corpus = None
    labels: Labels = None
    counts: dict = field(default_factory=dict)
    trace: Trace = Trace()
    raw: str = ""
    terms: dict = field(default_factory=dict)
    cached: bool = False
    max_date: str = ""


def run(arm, query, synonyms=None, workers=6, cache=False, runs_dir=None):
    """One pair through one arm.

    `synonyms` is only read by an expander told to replay a recorded list; an
    arm that expands live ignores it, which is the difference between a
    benchmark whose search must not move and a fresh question.

    `cache` reuses a saved run only when the arm's fingerprint matches, so a
    pipeline that has changed in any declared way is never answered from one
    that has not.

    The arm's max_date is set on the PubMed client here, for every search the
    run makes. The ceiling is process-wide, so concurrent runs must share it.
    """
    arm.validate()
    if cache:
        from . import store
        saved = store.load(query, arm, runs_dir)
        if saved is not None:
            saved.cached = True
            return saved
    from ..services import pubmed_service
    pubmed_service.MAX_PUBDATE = arm.max_date or None
    res = Result(query=query, arm=arm.name, fingerprint=arm.fingerprint(),
                 max_date=arm.max_date)
    clock = _stopwatch(res)

    terms = {}
    for block in arm.expand:
        with clock(block):
            terms.update(block.run(query, terms, synonyms))
    res.terms = terms

    # Passes that only look at the question go together; the merge order is the
    # order they are declared in, never the order they answer in, so an
    # article's rank records which query reached it first.
    eager = [s for s in arm.sources if not s.deferred]
    late = [s for s in arm.sources if s.deferred]
    corpus = Corpus()
    if eager:
        with ThreadPoolExecutor(max_workers=max(1, min(len(eager), workers))) as pool:
            found = list(pool.map(lambda s: s.run(query, terms), eager))
        for block, part in zip(eager, found):
            res.trace = res.trace.then(Step(block.kind, block.fingerprint(),
                                            note=f"{len(part)} articles"))
            corpus = corpus.merged(part)
    for block in late:
        if not block.applies(corpus, terms):
            res.trace = res.trace.then(Step(block.kind, block.fingerprint(),
                                            note="not applicable"))
            continue
        with clock(block):
            corpus = corpus.merged(block.run(query, terms))
    with clock(arm.dedupe):
        corpus = arm.dedupe.run(corpus)
    res.counts["found"] = len(corpus)

    chain = list(arm.annotate)
    chain.insert(min(arm.hydrate_at, len(chain)), arm.hydrate)
    labels = Labels()
    for block in chain:
        if block is arm.hydrate:
            with clock(block):
                # A hydrate block that only fetches abstracts needs the corpus
                # and nothing else. One that reads full text has to know which
                # gene it is looking for, and says so rather than having the
                # query pushed at every block that does not want it.
                corpus = (block.run(corpus, query)
                          if getattr(block, "wants_query", False)
                          else block.run(corpus))
            continue
        with clock(block):
            # Annotators see what the expansion produced, not just the query:
            # the grader that is told the disease's other names reads them here.
            part = block.run(query, corpus, terms)
        labels = labels.merged(part)
        gone = part.dropped()
        if gone:
            alive = corpus.keep(p for p in corpus.pmids if p not in gone)
            # An annotator may narrow the corpus but not empty it: a filter that
            # rejects everything has failed at its job, not answered it.
            corpus = alive if len(alive) else corpus
            res.counts[block.kind.lower()] = len(gone)

    with clock(arm.select):
        corpus = arm.select.run(corpus, labels)
    res.counts["read"] = len(corpus)
    res.corpus, res.labels = corpus, labels

    for attempt in range(1, READ_TRIES + 1):
        with clock(arm.read):
            text, _ = arm.read.run(query, corpus)
        if not arm.judge.garbled(text):
            break
    else:
        res.trace = res.trace.then(Step(arm.read.kind, arm.read.fingerprint(),
                                        note=f"answer garbled {READ_TRIES} times; "
                                             "the last one is read as well as it can be"))
    res.raw = text
    with clock(arm.judge):
        res.verdict = arm.judge.run(query, text, corpus)
    if cache:
        from . import store
        store.save(res, runs_dir)
    return res


def _stopwatch(res):
    """Times a block and files what it was, so the run can be read back."""
    class _Timer:
        def __init__(self, block):
            self.block = block

        def __enter__(self):
            from ..cost import capturing
            self.t = time.time()
            self.usage = capturing()
            self.usage.__enter__()
            return self

        def __exit__(self, *exc):
            self.usage.__exit__(*exc)
            spent = self.usage.summary(getattr(self.block, "model", None))
            res.trace = res.trace.then(Step(
                self.block.kind, self.block.fingerprint(),
                model=spent["model"] or getattr(self.block, "model", None),
                seconds=round(time.time() - self.t, 1),
                input_tokens=spent["input_tokens"],
                output_tokens=spent["output_tokens"],
                usd=spent["usd"],
                calls=spent["calls"]))
    return _Timer
