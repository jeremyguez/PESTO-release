"""Runs on disk, keyed by the question and by the pipeline that answered it.

A cached answer is only an answer to the same question if the same pipeline
produced it, and that is what an arm's fingerprint says. So it is in the
filename. A run made before the variant pass existed and a run made after are
two files, neither shadowing the other, and neither can be handed back for a
question it was not asked.

That replaces comparing a dozen recorded parameters one by one and hoping the
list was complete. It was not: the matcher that reads the Open Targets answer
was missing from it for months.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import fields

from ..config import PROJECT_ROOT

from .types import (Corpus, Distribution, FIELDS, Label, Labels, Query, Step,
                    Trace, Verdict)

_STEP_FIELDS = {f.name for f in fields(Step)}

RUNS = os.path.join(PROJECT_ROOT, "data", "runs", "flow")


def slug(query):
    text = f"{query.gene}_{query.phenotype}".lower()
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", text)).strip("_")[:110]


def path(query, arm, runs_dir=None):
    return os.path.join(runs_dir or RUNS,
                        f"{slug(query)}__{arm.fingerprint()}.json")


def save(result, runs_dir=None):
    from .run import Result  # noqa: F401  (documents what is being written)
    target = os.path.join(runs_dir or RUNS,
                          f"{slug(result.query)}__{result.fingerprint}.json")
    os.makedirs(os.path.dirname(target), exist_ok=True)
    payload = {
        "gene": result.query.gene, "phenotype": result.query.phenotype,
        "arm": result.arm, "fingerprint": result.fingerprint,
        "terms": {k: list(v) if isinstance(v, tuple) else v
                  for k, v in result.terms.items()},
        "counts": result.counts,
        "verdict": {"call": result.verdict.call,
                    "distribution": (result.verdict.distribution.as_dict()
                                     if result.verdict.distribution else None),
                    "justification": result.verdict.justification,
                    "pmids": list(result.verdict.pmids),
                    "detail": result.verdict.detail},
        "articles": result.corpus.as_dicts() if result.corpus else [],
        "labels": {p: {"relevance": lab.relevance, "kind": lab.kind,
                       "dropped": lab.dropped}
                   for p, lab in (result.labels.by_pmid.items()
                                  if result.labels else ())},
        "trace": result.trace.as_list(),
        "raw": result.raw,
    }
    with open(target, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return target


def load(query, arm, runs_dir=None):
    """The saved run for this exact pipeline, or nothing."""
    target = path(query, arm, runs_dir)
    if not os.path.exists(target):
        return None
    with open(target, encoding="utf-8") as fh:
        d = json.load(fh)
    from .run import Result
    v = d.get("verdict") or {}
    dist = v.get("distribution")
    return Result(
        query=Query(d["gene"], d["phenotype"]), arm=d["arm"],
        fingerprint=d["fingerprint"],
        verdict=Verdict(
            call=v.get("call", ""),
            distribution=(Distribution(
                established=dist.get("Established", 0),
                existing=dist.get("Existing", 0),
                hypothesized=dist.get("Hypothesized", 0),
                novel=dist.get("Novel", 0)) if dist else None),
            justification=v.get("justification", ""),
            pmids=tuple(v.get("pmids") or ()), detail=v.get("detail") or {}),
        corpus=Corpus.of(d.get("articles") or [], FIELDS),
        labels=Labels({p: Label(**lab) for p, lab in (d.get("labels") or {}).items()}),
        counts=d.get("counts") or {},
        trace=Trace(tuple(
            Step(**{k: v for k, v in s.items() if k in _STEP_FIELDS})
            for s in (d.get("trace") or []))),
        raw=d.get("raw", ""), terms=d.get("terms") or {})
