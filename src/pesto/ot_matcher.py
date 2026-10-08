"""Deciding whether the reported trait is already in the Open Targets answer.

This is the single comparison every Open Targets verdict rests on. A study
reports that a gene is associated with a trait; the database already links that
gene to a list of traits; this module decides whether the reported trait is in
that list.

The twenty traits nearest the reported one, by an embedding model, are sieved by
claude-haiku and then graded one by one by claude-opus on how each relates to
the reported trait. The model only says how each trait relates, having never
seen a score, and the arithmetic picks the deciding trait afterwards: the
highest-scoring trait graded the same trait earns its band, and failing that
the highest-scoring tie caps the pair at Hypothesized.
"""
import re

TAGGED_SHORTLIST = "tagged_shortlist"
METHODS = (TAGGED_SHORTLIST,)

# The model that grades the shortlist. A cheaper model is not a free choice
# here: read by claude-haiku, the same task scored 100% sensitivity at 52%
# specificity on the benchmark.
BROAD_MODEL = "claude-opus-5"
# The sieve in front of it. Its job is only to drop what is plainly unrelated,
# which is a judgement a cheap model makes as well as an expensive one; it is
# told repeatedly that being unsure is a reason to keep.
GATE_MODEL = "claude-haiku-4-5"


# How each trait may relate to the reported one. Only the closer grades
# count as a tie and cap at Hypothesized: a first-degree parent or child,
# a genetic correlate that would often share a signal, or a specific or
# moderate symptom. Distant parents, vague resemblance and nonspecific
# markers are ignored, same as DIFFERENT, so the pair stays Novel.
TAGS = ("EXACT",
        "PARENT_FIRST", "PARENT_DISTANT",
        "CHILD_FIRST", "CHILD_DISTANT",
        "CLOSE_VERY", "CLOSE_MODERATE", "CLOSE_VAGUE",
        "SYMPTOM_SPECIFIC", "SYMPTOM_MODERATE", "SYMPTOM_NONSPECIFIC",
        "DIFFERENT")
TIES = (
    "PARENT_FIRST", "CHILD_FIRST",
    "CLOSE_VERY", "CLOSE_MODERATE",
    "SYMPTOM_SPECIFIC", "SYMPTOM_MODERATE",
)


def resolve_method(method=None):
    """The matcher to use. There is one."""
    m = method or TAGGED_SHORTLIST
    if m not in METHODS:
        raise ValueError(f"unknown OT matcher {m!r}, expected one of {METHODS}")
    return m


def _render(traits):
    """The shortlist as the reader sees it: numbered, and without the scores.

    Withholding the scores is the point of the design. The reader is asked how
    two traits relate, which is a question about biology; showing it that the
    database scores one of them 0.8 invites it to answer a different question,
    about how confident the database is.
    """
    return "\n".join(f"{i}. {t.get('name', '')}"
                     for i, t in enumerate(traits, 1))


def _parse_tags(text, n):
    """One tag per entry, by number. An unreadable line leaves its entry None.

    Entries the reader never spoke about are left unset rather than guessed at,
    and they simply do not count as ties: a missing tag can only weaken the
    verdict, never strengthen it.
    """
    tags, reasons = [None] * n, [""] * n
    for line in (text or "").splitlines():
        m = re.match(r"\s*(\d+)\s*[.|)]?\s*\|?\s*([A-Z_]+)\s*(?:\|\s*(.*))?$",
                     line.strip())
        if not m:
            continue
        i, tag = int(m.group(1)) - 1, m.group(2).upper()
        if 0 <= i < n and tag in TAGS:
            tags[i] = tag
            reasons[i] = (m.group(3) or "").strip()
    return tags, reasons


def _parse_gate(text, n):
    """The entries the sieve wants removed, as zero-based indices.

    Fail-safe on purpose: anything that does not match the expected line leaves
    the shortlist whole. Dropping a trait that mattered costs a verdict, and
    keeping a useless one costs one line of a prompt.
    """
    m = re.search(r"DROP:\s*(.*)", text or "", re.I)
    if not m or m.group(1).strip().upper().startswith("NONE"):
        return set()
    return {i for i in (int(d) - 1 for d in re.findall(r"\d+", m.group(1)))
            if 0 <= i < n}


def _tagged_verdict(traits, tags):
    """Which trait decides, and whether its band may be earned in full.

    Two channels, because the two kinds of hit warrant different claims. A trait
    the reader called the same trait carries whatever band its score earns. A
    trait that is merely tied — a parent, a component, a symptom, a correlate —
    caps at Hypothesized however high the database scores it, because the score
    measures the evidence for *that* association and not for the one asked
    about.

    Within each channel the highest-scoring qualifying trait wins, rather than
    whichever trait the model would have chosen to talk about.
    """
    exact = [t for t, g in zip(traits, tags) if g == "EXACT"]
    if exact:
        best = max(exact, key=lambda t: t.get("score") or 0.0)
        return best, "exact", None
    ties = [t for t, g in zip(traits, tags) if g in TIES]
    if ties:
        best = max(ties, key=lambda t: t.get("score") or 0.0)
        return best, "tie", "Hypothesized"
    return None, "none", None


def _slim_trait(t):
    return {"disease_id": t.get("disease_id"),
            "name": t.get("name", ""),
            "score": t.get("score")}


def replay_tagged(record, apply_gate=None):
    """Recompute the deciding trait from a saved agent record and current TIES.

    The agents are not called. Changing TIES, the score bands or whether the
    gate is applied is therefore free, provided the record graded the full
    shortlist (so dropped traits still have tags).
    """
    from . import config
    apply_gate = config.OT_GATE if apply_gate is None else apply_gate
    short = [dict(t) for t in (record.get("shortlist") or [])]
    graded = list(record.get("graded") or [])
    by_id = {g.get("disease_id"): g for g in graded}
    dropped = set(record.get("gate_dropped_ids") or [])
    used = [t for t in short
            if not (apply_gate and t.get("disease_id") in dropped)]
    # Same fail-safe as the live sieve: dropping everything is a misread.
    if apply_gate and dropped and not used:
        used = short
    tags = [(by_id.get(t.get("disease_id")) or {}).get("tag") for t in used]
    hit, channel, cap = _tagged_verdict(used, tags)
    deciding_tag, why = None, ""
    if hit is not None:
        g = by_id.get(hit.get("disease_id")) or {}
        deciding_tag = g.get("tag")
        why = g.get("why") or ""
    if hit is None:
        reason = (f"none of the {len(used)} nearest traits is the reported "
                  f"trait or is tied to it")
    else:
        stated = why or (deciding_tag or "").lower().replace("_", " ")
        reason = (f"{hit.get('name', '')} is {deciding_tag} ({stated})"
                  if channel == "exact" else
                  f"no trait is the reported trait; the closest tie is "
                  f"{hit.get('name', '')}, {deciding_tag} ({stated}), which caps "
                  f"the verdict at Hypothesized")
    return {
        "method": TAGGED_SHORTLIST,
        "matched_id": hit.get("disease_id") if hit else None,
        "basis": hit.get("name", "") if hit else "",
        "deciding_tag": deciding_tag,
        "reason": reason,
        "channel": channel,
        "cap": cap,
        "read": len(used),
        "of": record.get("of"),
        "gate_dropped": len(dropped) if apply_gate else 0,
        "tags": {g.get("name", ""): g.get("tag")
                 for g in graded if g.get("tag")},
    }


def _match_tagged(gene, phenotype, traits, model=None, temperature=0.0,
                  encoder=None, top_k=None, gate=None, reasons=None):
    """Shortlist, sieve, then grade every shortlisted trait one by one.

    The sieve is recorded, not applied before the reader: every nearest trait
    is tagged so a later change to TIES or to whether the gate counts can be
    replayed from disk without another model call.
    """
    from . import config, ot_shortlist
    from .services.llm_service import call_llm_with_usage
    from .utils.helpers import load_prompt

    gate = config.OT_GATE if gate is None else gate
    reasons = config.OT_REASONS if reasons is None else reasons

    short = ot_shortlist.shortlist(phenotype, traits, k=top_k, name=encoder)
    # Two calls to two models at prices six times apart, so the budget guard is
    # given their sum while the audit keeps them separate. Adding them and then
    # trying to price the total would be a guess about which model spent what.
    spent = {"input_tokens": 0, "output_tokens": 0}
    tokens = {}

    def spend(res, which):
        spent["input_tokens"] += res.get("input_tokens") or 0
        spent["output_tokens"] += res.get("output_tokens") or 0
        tokens[which + "_in"] = res.get("input_tokens") or 0
        tokens[which + "_out"] = res.get("output_tokens") or 0
        return res

    # Room for one line per trait on top of the model's thinking. 4096 fits the
    # twenty nearest; read whole (--ot-encoder none), a gene's answer runs to
    # several hundred traits, and a cut-off answer carries no tag at all.
    per_trait = 128 if reasons else 64
    room = min(64000, max(4096, 2048 + per_trait * len(short)))

    gate_text = ""
    dropped_ids = []
    if gate and short:
        template = load_prompt("ot_tagged_gate_prompt")
        if not template:
            raise RuntimeError("ot_tagged_gate_prompt.txt not found")
        res = spend(call_llm_with_usage(
            template.format(gene=gene, phenotype=phenotype,
                            traits=_render(short)),
            GATE_MODEL, temperature, agent_name="ot_tagged_gate_agent",
            max_tokens=room), "gate")
        gate_text = res.get("text") or ""
        cut = _parse_gate(gate_text, len(short))
        # Never everything: a sieve that empties the list has misread the task,
        # and an empty list would be read downstream as a novel association.
        if cut and len(cut) < len(short):
            dropped_ids = [short[i].get("disease_id") for i in sorted(cut)]

    name = config.ot_match_prompt_name(reasons=reasons)
    template = load_prompt(name)
    if not template:
        raise RuntimeError(f"{name}.txt not found")
    res = spend(call_llm_with_usage(
        template.format(gene=gene, phenotype=phenotype,
                        traits=_render(short)),
        model or BROAD_MODEL, temperature,
        agent_name="ot_tagged_match_agent", max_tokens=room), "read")
    tag_text = res.get("text") or ""
    tags, why = _parse_tags(tag_text, len(short))
    # A cut-off or empty answer is a failed reading, not an answer: read as
    # one, the traits it never reached would count as unrelated and pull the
    # pair towards Novel.
    if short and (res.get("stop_reason") == "max_tokens" or not any(tags)):
        raise RuntimeError(
            f"the Open Targets grader tagged {sum(1 for t in tags if t)} of "
            f"{len(short)} traits (stop reason: {res.get('stop_reason')})")
    graded = [{"disease_id": t.get("disease_id"),
               "name": t.get("name", ""),
               "score": t.get("score"),
               "tag": tags[i],
               "why": why[i]}
              for i, t in enumerate(short)]
    record = {
        "shortlist": [_slim_trait(t) for t in short],
        "gate_text": gate_text,
        "gate_dropped_ids": dropped_ids,
        "tag_text": tag_text,
        "graded": graded,
        "of": len(traits),
    }
    derived = replay_tagged(record, apply_gate=gate)
    derived.update(record)
    derived["tokens"] = tokens
    derived["usage"] = dict(spent, text=tag_text)
    return derived


def match(gene, phenotype, traits, method=None, model=None, temperature=0.0,
          encoder=None, top_k=None, gate=None, reasons=None):
    """Which trait in the answer, if any, records the reported association.

    `traits` are the entries the database returned, each a dict with
    `disease_id`, `name` and `score`, already past the score floor.

    Returns a dict carrying `matched_id`, the trait that decides the verdict or
    None, together with the basis and the reasoning, so that the verdict can be
    audited without re-running the matcher. The tagged matcher adds `cap`, the
    strongest verdict its hit may earn, which is the one field a caller must not
    ignore: without it a tie would be read as the association itself.
    """
    chosen = resolve_method(method)
    if not traits:
        return {"method": chosen, "matched_id": None, "basis": "", "cap": None,
                "reason": "the database returned no trait above the score floor",
                "usage": None}
    return _match_tagged(gene, phenotype, traits, model=model,
                         temperature=temperature, encoder=encoder,
                         top_k=top_k, gate=gate, reasons=reasons)
