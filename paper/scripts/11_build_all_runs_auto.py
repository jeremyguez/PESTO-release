#!/usr/bin/env python3
"""Build the figure-1 table from the auto current / current-general cohort run.

Literature verdicts come from results/auto_cohort/. Burden statistics come
from the source tables. Open Targets is joined from the tight rematch
(results/arm_auto_aou_brava_ot_tight.tsv) when it exists, else from the
previous all_runs.tsv. The plotted `verdict` is lit_mean, as in figure 2a–c.

  python3 scripts/11_build_all_runs_auto.py
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pesto.cli import mean_call  # noqa: E402
from pesto.flow.store import slug  # noqa: E402
from pesto.flow.types import Query  # noqa: E402

from importlib import import_module  # noqa: E402
_aggregate = import_module("00_aggregate_runs")

RUNS = os.path.join(ROOT, "results", "auto_cohort")
OLD = os.path.join(ROOT, "results", "all_runs.tsv")
OT_TIGHT = os.path.join(ROOT, "results", "arm_auto_aou_brava_ot_tight.tsv")
OUT = os.path.join(ROOT, "results", "all_runs_auto.tsv")

# AoU and BRaVa name some phenotypes differently ("BMI" and "Body mass index",
# "high density lipoprotein" and "HDL cholesterol"): the same pair, run twice
# under two names, would otherwise be counted twice.
ALIASES = os.path.join(ROOT, "data", "phenotype_aliases.tsv")
with open(ALIASES, encoding="utf-8") as _fh:
    PHENOTYPE_ALIASES = {r["alias"]: r["canonical"]
                         for r in csv.DictReader(_fh, delimiter="\t")}


def pair_key(gene, phenotype):
    ph = re.sub(r"[^a-z0-9]+", " ", phenotype.lower()).strip()
    return gene.upper(), PHENOTYPE_ALIASES.get(ph, ph)


def drop_shared_brava(rows):
    """Keep one row per pair significant in both cohorts, on the AoU side,
    as classify() already does for pairs whose names match exactly."""
    aou = {pair_key(r["gene"], r["phenotype"]) for r in rows if r["source"] == "AoU"}
    kept, dropped = [], []
    for r in rows:
        if r["source"] == "BRAVA" and pair_key(r["gene"], r["phenotype"]) in aou:
            dropped.append(r)
        else:
            kept.append(r)
    for r in dropped:
        print(f"dropped BRAVA duplicate of an AoU pair: {r['gene']} / {r['phenotype']}")
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ot-table", default=OT_TIGHT,
                    help="Open Targets rematch summary to join (default: tight).")
    ap.add_argument("--out", default=OUT,
                    help="output table (default: results/all_runs_auto.tsv).")
    args = ap.parse_args()

    aou, brava = _aggregate.load_sources()
    old = {}
    with open(OLD, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            old[(r["gene"], r["phenotype"])] = r
    tight = {}
    if os.path.exists(args.ot_table):
        with open(args.ot_table, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tight[(r["gene"], r["phenotype"])] = r

    rows = []
    for path in sorted(glob.glob(os.path.join(RUNS, "*__*.json"))):
        d = json.load(open(path, encoding="utf-8"))
        gene, phenotype = d["gene"], d["phenotype"]
        dist = (d.get("verdict") or {}).get("distribution") or {}
        source, stats = _aggregate.classify(gene, phenotype, aou, brava)
        prev = old.get((gene, phenotype), {})
        ot = tight.get((gene, phenotype), {})
        argmax = (d.get("verdict") or {}).get("call") or ""
        lit = mean_call(dist) or argmax
        rows.append({
            "gene": gene,
            "phenotype": phenotype,
            "source": source,
            "arm": d.get("arm", ""),
            "fingerprint": d.get("fingerprint", ""),
            "model": "auto-opus",
            "verdict": lit,
            "lit_argmax": argmax,
            "lit_mean": lit,
            "p_established": dist.get("Established", ""),
            "p_existing": dist.get("Existing", ""),
            "p_hypothesized": dist.get("Hypothesized", ""),
            "p_novel": dist.get("Novel", ""),
            "ot_verdict": (ot.get("open_targets_verdict")
                           or prev.get("ot_verdict", "")),
            "ot_max_score": (ot.get("open_targets_max_score")
                             if ot else prev.get("ot_max_score", "")),
            "ot_matched_traits": ot.get("ot_trait") or prev.get("ot_matched_traits", ""),
            "ot_matcher": (ot.get("ot_matcher")
                           or prev.get("ot_matcher", "")),
            "ot_channel": ot.get("ot_channel", ""),
            "ot_tag": ot.get("ot_tag", ""),
            "n_pubmed": (d.get("counts") or {}).get("found", ""),
            "n_read": (d.get("counts") or {}).get("read", ""),
            "n_cited": len((d.get("verdict") or {}).get("pmids") or []),
            "cited_pmids": ",".join((d.get("verdict") or {}).get("pmids") or []),
            "burden_pvalue": stats.get("pvalue", ""),
            "burden_beta": stats.get("beta", ""),
            "n_cases": stats.get("n_cases", ""),
            "n_controls": stats.get("n_controls", ""),
            "justification": ((d.get("verdict") or {}).get("justification") or ""
                              ).replace("\t", " ").replace("\n", " "),
            "run_file": os.path.basename(path),
        })

    rows = drop_shared_brava(rows)
    rows.sort(key=lambda r: (r["source"], r["gene"], r["phenotype"]))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(f"wrote {len(rows)} pairs to {os.path.relpath(args.out, ROOT)}")
    print("source", Counter(r["source"] for r in rows))
    print("arm   ", Counter(r["arm"] for r in rows))
    print("verdict")
    for src in ("AoU", "BRAVA"):
        sub = [r for r in rows if r["source"] == src]
        print(f"  {src}", Counter(r["verdict"] for r in sub))


if __name__ == "__main__":
    main()
