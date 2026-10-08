#!/usr/bin/env python3
"""Re-grade saved Open Targets shortlists under a different match prompt.

Traits, BioLORD ranking and the Haiku sieve stay frozen. Only the Opus
tagging call is new. Sidecars are written to a separate directory so the
source run is not overwritten.

  OT_MATCH_PROMPT=v5 python3 scripts/17_replay_ot_prompt.py \
      --src results/auto_cohort_ot_tight \
      --dst results/auto_cohort_ot_v5 \
      --out results/arm_auto_aou_brava_ot_v5.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))

from pesto import config  # noqa: E402
from pesto.ot_matcher import (  # noqa: E402
    BROAD_MODEL, _parse_tags, _render, replay_tagged,
)
from pesto.scoring import open_targets_verdict  # noqa: E402
from pesto.services.llm_service import call_llm_with_usage  # noqa: E402
from pesto.utils.helpers import load_prompt  # noqa: E402

SRC = os.path.join(ROOT, "results", "auto_cohort_ot_tight")
DST = os.path.join(ROOT, "results", "auto_cohort_ot_v5")
OUT = os.path.join(ROOT, "results", "arm_auto_aou_brava_ot_v5.tsv")


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name or "") or "gene"


def sidecar_name(gene, phenotype):
    return f"open_targets_match_{safe(gene)}_{safe(phenotype)}.json"


def decide(rec):
    derived = replay_tagged(rec, apply_gate=True)
    traits = []
    for t in rec.get("traits") or []:
        row = dict(t)
        row["matched"] = bool(derived.get("matched_id")
                              and row.get("disease_id") == derived["matched_id"])
        traits.append(row)
    verdict, score = open_targets_verdict({"traits": traits, "match": derived})
    return derived, verdict, score


def already_replayed(path, src_tag):
    if not os.path.exists(path):
        return False
    try:
        rec = json.load(open(path, encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return rec.get("replayed_from") == src_tag and bool(rec.get("graded"))


def one(src_path, dst_dir, template, prompt, src_tag):
    rec = json.load(open(src_path, encoding="utf-8"))
    gene = rec.get("gene") or ""
    phenotype = rec.get("phenotype") or ""
    dest = os.path.join(dst_dir, sidecar_name(gene, phenotype))
    if already_replayed(dest, src_tag):
        derived, verdict, score = decide(json.load(open(dest, encoding="utf-8")))
        return row_from(gene, phenotype, derived, verdict, score, cached=True)

    short = rec.get("shortlist") or []
    res = call_llm_with_usage(
        template.format(gene=gene, phenotype=phenotype, traits=_render(short)),
        BROAD_MODEL, 0.0, agent_name="ot_tagged_match_agent")
    text = res.get("text") or ""
    tags, whys = _parse_tags(text, len(short))
    graded = [{"disease_id": t.get("disease_id"),
               "name": t.get("name", ""),
               "score": t.get("score"),
               "tag": tags[i],
               "why": whys[i]}
              for i, t in enumerate(short)]
    agent = config.ot_agent_version(prompt=prompt)
    fresh = {
        "gene": gene,
        "phenotype": phenotype,
        "ensembl_id": rec.get("ensembl_id"),
        "agent_version": agent,
        "ot_match_version": agent + (":gate" if config.OT_GATE else ""),
        "traits": rec.get("traits") or [],
        "shortlist": short,
        "gate_text": rec.get("gate_text") or "",
        "gate_dropped_ids": rec.get("gate_dropped_ids") or [],
        "tag_text": text,
        "graded": graded,
        "of": rec.get("of"),
        "tokens": {
            "read_in": res.get("input_tokens") or 0,
            "read_out": res.get("output_tokens") or 0,
        },
        "articles": rec.get("articles") or [],
        "replayed_from": src_tag,
    }
    os.makedirs(dst_dir, exist_ok=True)
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(fresh, fh, indent=2, default=str)
        fh.write("\n")
    derived, verdict, score = decide(fresh)
    return row_from(gene, phenotype, derived, verdict, score, cached=False)


def row_from(gene, phenotype, derived, verdict, score, cached):
    return {
        "gene": gene,
        "phenotype": phenotype,
        "open_targets_verdict": verdict,
        "open_targets_max_score": "" if score is None else score,
        "ot_channel": derived.get("channel") or "",
        "ot_tag": derived.get("deciding_tag") or "",
        "ot_trait": derived.get("basis") or "",
        "ot_cap": derived.get("cap") or "",
        "ot_matcher": derived.get("method") or "",
        "cached": cached,
        "error": "",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=SRC)
    ap.add_argument("--dst", default=DST)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--prompt", default=os.environ.get("OT_MATCH_PROMPT", "v5"))
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    prompt_file = config.ot_match_prompt_name(prompt=args.prompt)
    template = load_prompt(prompt_file)
    if not template:
        raise SystemExit(f"pesto: {prompt_file}.txt not found")

    src_tag = os.path.basename(os.path.abspath(args.src))
    paths = sorted(
        os.path.join(args.src, n) for n in os.listdir(args.src)
        if n.startswith("open_targets_match_") and n.endswith(".json"))
    os.makedirs(args.dst, exist_ok=True)
    print(f"{len(paths)} records, {args.workers} workers, "
          f"prompt={args.prompt} ({prompt_file}), "
          f"{src_tag} -> {os.path.basename(args.dst)}", flush=True)

    rows, done, called = [], 0, 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(one, p, args.dst, template, args.prompt, src_tag): p
                for p in paths}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            done += 1
            if not row["cached"]:
                called += 1
            mark = "cache" if row["cached"] else "opus"
            print(f"{done:4d}/{len(paths)}  {row['gene']:12s}  "
                  f"{row['open_targets_verdict']:13s}  "
                  f"{row['ot_tag'] or '-':18s}  {mark}",
                  flush=True)

    rows.sort(key=lambda r: (r["gene"], r["phenotype"]))
    fields = ["gene", "phenotype", "open_targets_verdict",
              "open_targets_max_score", "ot_channel", "ot_tag",
              "ot_trait", "ot_cap", "ot_matcher", "error"]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}  ({called} Opus calls, "
          f"{done - called} reused)")
    print("OT", dict(Counter(r["open_targets_verdict"] for r in rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
