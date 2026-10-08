"""The Open Targets verdict: a matched trait's score read off by threshold."""
from __future__ import annotations

VERDICTS = ["Novel", "Hypothesized", "Existing", "Established"]


# The database has already done its weighing and published the result. All
# that is left is to read its number off the one trait that records the
# association, which is why this is a threshold and not a model. The thresholds
# are the database's own quartiles of what it calls a strong association, and
# they have not moved since the first benchmark.
OT_HYPOTHESIZED_BELOW = 0.2
OT_EXISTING_UPTO = 0.5


def open_targets_band(score):
    """What a matched trait's score is worth, before any cap."""
    if score is None:
        return "Novel"
    if score < OT_HYPOTHESIZED_BELOW:
        return "Hypothesized"
    if score <= OT_EXISTING_UPTO:
        return "Existing"
    return "Established"


def open_targets_verdict(open_targets_traits):
    """The database's verdict on the pair, and the score it rests on.

    Returns (verdict, score). "Error" when the lookup itself failed, which is
    not the same as the database having nothing: a network fault must never be
    read as a novel association.

    The matcher names at most one trait, so the loop below finds one score. It
    is written as a maximum anyway because a saved run from before that rule
    can hold several, and rescoring an old run must not depend on dictionary
    order.

    A `cap` in the match detail is honoured. The tagged matcher sets it when the
    trait it found is only tied to the reported one — a parent, a component, a
    symptom, a correlate. The database's score then measures the evidence for
    that neighbouring association and not for this one, so it may suggest the
    association and no more, whatever its size. Absent from older runs, where
    it correctly does nothing.
    """
    if not open_targets_traits:
        return "Novel", 0.0
    if open_targets_traits.get("error"):
        return "Error", 0.0

    best = None
    for t in open_targets_traits.get("traits") or []:
        if not t.get("matched"):
            continue
        s = t.get("score")
        if s is not None and (best is None or s > best):
            best = s
    if best is None:
        return "Novel", 0.0

    verdict = open_targets_band(best)
    cap = (open_targets_traits.get("match") or {}).get("cap")
    if cap in VERDICTS and VERDICTS.index(cap) < VERDICTS.index(verdict):
        verdict = cap
    return verdict, best
