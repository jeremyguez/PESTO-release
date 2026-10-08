#!/usr/bin/env python3
"""Re-match Open Targets on the AoU / BRAVA pairs. Literature is not touched.

Writes match sidecars under results/auto_cohort_ot_tight/ and a summary TSV.
The default tagged prompt (tight) and current TIES apply.

  python3 scripts/16_rematch_ot_cohort.py --workers 20
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))

from pesto.open_targets import fetch_and_save_open_targets_traits  # noqa: E402
from pesto.scoring import open_targets_verdict  # noqa: E402

BENCH = os.path.join(ROOT, "results", "bench_aou_brava_auto.tsv")
RUNS = os.path.join(ROOT, "results", "auto_cohort_ot_tight")
OUT = os.path.join(ROOT, "results", "arm_auto_aou_brava_ot_tight.tsv")


def one(gene, phenotype, runs):
    chatter = io.StringIO()
    try:
        with contextlib.redirect_stdout(chatter):
            traits = fetch_and_save_open_targets_traits(
                gene, runs, phenotype, [])
        match = (traits or {}).get("match") or {}
        verdict, score = open_targets_verdict(traits)
        return {
            "gene": gene,
            "phenotype": phenotype,
            "open_targets_verdict": verdict,
            "open_targets_max_score": "" if score is None else score,
            "ot_channel": match.get("channel") or "",
            "ot_tag": match.get("deciding_tag") or "",
            "ot_trait": match.get("basis") or "",
            "ot_cap": match.get("cap") or "",
            "ot_matcher": match.get("method") or "",
            "error": (traits or {}).get("error") or "",
        }
    except Exception as exc:
        return {
            "gene": gene, "phenotype": phenotype,
            "open_targets_verdict": "Error",
            "open_targets_max_score": "",
            "ot_channel": "", "ot_tag": "", "ot_trait": "", "ot_cap": "",
            "ot_matcher": "", "error": f"{type(exc).__name__}: {exc}",
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--bench", default=BENCH)
    ap.add_argument("--runs-dir", default=RUNS)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    with open(args.bench, encoding="utf-8") as fh:
        pairs = [(r["gene"], r["phenotype"]) for r in csv.DictReader(fh, delimiter="\t")]
    os.makedirs(args.runs_dir, exist_ok=True)
    print(f"{len(pairs)} pairs, {args.workers} workers, OT only", flush=True)

    # Cap torch to one core per worker before any encode runs. Torch defaults to
    # half the cores per operation, so twenty workers each launching a ten-thread
    # encode oversubscribe a twenty-core box tenfold (load average 200), which
    # freezes the machine while the real per-pair cost is the Opus call waiting
    # on the network, not the encode. One core a worker keeps the box responsive
    # and costs the encode nothing: it embeds a few hundred short trait names.
    try:
        import torch  # noqa: E402
        torch.set_num_threads(max(1, (os.cpu_count() or 1) // max(1, args.workers)))
    except ImportError:
        pass

    # Warm the shortlist encoder once, before the pool opens, so the workers do
    # not all miss the cache together and each load their own 439 MB copy.
    from pesto import ot_shortlist  # noqa: E402
    if ot_shortlist.encoder() is not None:
        print("encoder warmed", flush=True)

    rows, done = [], 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(one, g, p, args.runs_dir): (g, p) for g, p in pairs}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            done += 1
            print(f"{done:4d}/{len(pairs)}  {row['gene']:12s}  "
                  f"{row['open_targets_verdict']:13s}  "
                  f"{row['ot_tag'] or '-'}",
                  flush=True)

    order = {(g, p): i for i, (g, p) in enumerate(pairs)}
    rows.sort(key=lambda r: order[(r["gene"], r["phenotype"])])
    fields = ["gene", "phenotype", "open_targets_verdict",
              "open_targets_max_score", "ot_channel", "ot_tag",
              "ot_trait", "ot_cap", "ot_matcher", "error"]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(f"wrote {args.out}")
    print("OT", dict(Counter(r["open_targets_verdict"] for r in rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
