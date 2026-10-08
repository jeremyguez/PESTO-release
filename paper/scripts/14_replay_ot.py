#!/usr/bin/env python3
"""Re-apply current Open Targets rules to saved agent records.

No model is called. Each `open_targets_match_*.json` already holds the
shortlist, the gate text, every tag and the raw agent replies. This
script walks those files, runs the live TIES and score bands, and writes
the verdicts back onto a pesto TSV plus a long tag table.

  python3 scripts/14_replay_ot.py \
    --runs-dir results/auto_fig1c_absent50_hgnc \
    --pesto results/bench_fig1c_absent50_hgnc_pesto.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))

from pesto import ot_matcher  # noqa: E402
from pesto.scoring import open_targets_verdict  # noqa: E402


def load_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def replay_one(record, apply_gate=None):
    derived = ot_matcher.replay_tagged(record, apply_gate=apply_gate)
    traits = []
    for t in record.get("traits") or []:
        row = dict(t)
        row["matched"] = bool(derived.get("matched_id")
                              and row.get("disease_id") == derived["matched_id"])
        traits.append(row)
    verdict, score = open_targets_verdict({"traits": traits, "match": derived})
    return derived, verdict, score


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--pesto", required=True)
    p.add_argument("--no-gate", action="store_true",
                   help="ignore the saved sieve (needs tags on the full shortlist)")
    args = p.parse_args()

    apply_gate = False if args.no_gate else None
    records = []
    for name in sorted(os.listdir(args.runs_dir)):
        if not name.startswith("open_targets_match_") or not name.endswith(".json"):
            continue
        rec = load_json(os.path.join(args.runs_dir, name))
        if not rec.get("graded"):
            continue
        records.append(rec)

    by_pair = {}
    tag_rows = []
    for rec in records:
        derived, verdict, score = replay_one(rec, apply_gate=apply_gate)
        key = (rec.get("gene", ""), rec.get("phenotype", ""))
        by_pair[key] = {
            "open_targets_verdict": verdict,
            "open_targets_max_score": "" if score is None else score,
            "ot_channel": derived.get("channel") or "",
            "ot_tag": derived.get("deciding_tag") or "",
            "ot_trait": derived.get("basis") or "",
            "ot_cap": derived.get("cap") or "",
        }
        dropped = set(rec.get("gate_dropped_ids") or [])
        for g in rec.get("graded") or []:
            tag_rows.append({
                "gene": rec.get("gene", ""),
                "phenotype": rec.get("phenotype", ""),
                "disease_id": g.get("disease_id") or "",
                "trait": g.get("name") or "",
                "score": "" if g.get("score") is None else g.get("score"),
                "tag": g.get("tag") or "",
                "why": g.get("why") or "",
                "gated": g.get("disease_id") in dropped,
                "deciding": g.get("disease_id") == derived.get("matched_id"),
            })

    if os.path.exists(args.pesto):
        with open(args.pesto, encoding="utf-8") as fh:
            pesto = list(csv.DictReader(fh, delimiter="\t"))
        extra = ["ot_channel", "ot_tag", "ot_trait", "ot_cap"]
        fields = list(pesto[0].keys()) if pesto else []
        for col in extra:
            if col not in fields:
                fields.append(col)
        n = 0
        for row in pesto:
            hit = by_pair.get((row.get("gene", ""), row.get("phenotype", "")))
            if not hit:
                continue
            row.update(hit)
            n += 1
        with open(args.pesto, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(pesto)
        print(f"updated {n} pesto rows in {args.pesto}")
    else:
        print(f"pesto table missing: {args.pesto}", file=sys.stderr)

    tags_out = os.path.splitext(args.pesto)[0] + "_ot_tags.tsv"
    tag_fields = ["gene", "phenotype", "disease_id", "trait", "score",
                  "tag", "why", "gated", "deciding"]
    with open(tags_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=tag_fields, delimiter="\t")
        w.writeheader()
        w.writerows(tag_rows)
    print(f"wrote {len(tag_rows)} tags for {len(records)} pairs to {tags_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
