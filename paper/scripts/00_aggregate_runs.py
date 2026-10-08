#!/usr/bin/env python3
"""Flatten the cached novelty runs into one analysis-ready table.

Reads every JSON under data/runs/novel/, labels each (gene, phenotype) pair with
its source cohort, and joins the burden statistics from the AoU and BRAVA result
tables. Writes results/all_runs.tsv.

The atrial-fibrillation re-run of May 2026 is not part of the cache: those pairs
match neither source table and their provenance is unknown (see QUESTIONS.md).
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from datetime import datetime

import pandas as pd

os.environ.setdefault("PESTO_PROJECT_ROOT",
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pesto.config import AOU_FILE, BRAVA_FILE, NOVEL_RUNS_DIR, RESULTS_DIR  # noqa: E402


def load_sources():
    aou, brava = {}, {}
    with open(AOU_FILE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            aou[(r["gene_symbol"], r["description"])] = {
                "pvalue": r.get("META_Pvalue_Burden", ""),
                "beta": r.get("META_BETA_Burden", ""),
                "n_cases": r.get("n_cases", ""),
                "n_controls": r.get("n_controls", ""),
                "annotation": r.get("annotation", ""),
            }
    with open(BRAVA_FILE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            info = {"pvalue": r.get("min_pvalue", ""), "beta": "",
                    "n_cases": "", "n_controls": "", "annotation": ""}
            brava[(r["external_gene_name"], r["phenotype_full"])] = info
            brava[(r["external_gene_name"], r["phenotype"])] = info
    return aou, brava


def classify(gene, phenotype, aou, brava):
    if "autism" in (phenotype or "").lower():
        return "ASD", {}
    if (gene, phenotype) in aou:
        return "AoU", aou[(gene, phenotype)]
    if (gene, phenotype) in brava:
        return "BRAVA", brava[(gene, phenotype)]
    return "OTHER", {}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs-dir", default=NOVEL_RUNS_DIR)
    ap.add_argument("--out", default=os.path.join(RESULTS_DIR, "all_runs.tsv"))
    args = ap.parse_args()

    aou, brava = load_sources()
    rows, skipped = [], 0

    for path in sorted(glob.glob(os.path.join(args.runs_dir, "*_novel_run*.json"))):
        try:
            d = json.load(open(path))
        except Exception:
            skipped += 1
            continue

        gene = d.get("gene_name")
        phenotype = d.get("phenotype")
        if not gene or not phenotype:
            skipped += 1
            continue

        final = d.get("final_output") or {}
        api = d.get("api_response") or {}
        params = d.get("run_parameters") or {}
        payload = api if api else final

        ot = payload.get("open_targets_traits")
        ot_matched = 0
        if isinstance(ot, dict):
            ot_matched = sum(1 for t in (ot.get("traits") or []) if t.get("matched"))

        # Runs rematched by scripts/40 carry the trait that decides the verdict
        # and the matcher that named it; runs still on the original substring
        # rule carry neither, which is how the two are told apart downstream.
        match = payload.get("open_targets_match") or {}

        source, stats = classify(gene, phenotype, aou, brava)
        cited = payload.get("cited_articles") or final.get("cited_articles") or []

        rows.append({
            "gene": gene,
            "phenotype": phenotype,
            "source": source,
            "model": params.get("model", ""),
            "verdict": payload.get("verdict") or final.get("verdict") or "",
            "ot_verdict": payload.get("open_targets_verdict", ""),
            "ot_max_score": payload.get("open_targets_max_score", ""),
            "ot_matched_traits": ot_matched,
            "ot_matcher": match.get("method", "substring"),
            "ot_basis": match.get("basis", ""),
            "ot_verdict_substring": payload.get(
                "open_targets_verdict_substring", ""),
            "n_pubmed": len(payload.get("pubmed_pmids") or d.get("pubmed_pmids") or []),
            "n_cited": len(cited),
            "cited_pmids": ", ".join(
                str(a["pmid"]) if isinstance(a, dict) else str(a) for a in cited),
            "burden_pvalue": stats.get("pvalue", ""),
            "burden_beta": stats.get("beta", ""),
            "n_cases": stats.get("n_cases", ""),
            "n_controls": stats.get("n_controls", ""),
            "justification": (payload.get("justification") or "").replace("\t", " ").replace("\n", " "),
            "run_date": datetime.fromtimestamp(os.path.getmtime(path)).date().isoformat(),
            "run_file": os.path.basename(path),
        })

    df = pd.DataFrame(rows)
    # One row per pair: keep the most recent run when a pair was assessed twice.
    df = df.sort_values("run_date").drop_duplicates(["gene", "phenotype"], keep="last")
    df = df.sort_values(["source", "gene", "phenotype"]).reset_index(drop=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote {len(df)} pairs to {args.out} ({skipped} unreadable files skipped)")
    print("\nBy source:")
    print(df["source"].value_counts().to_string())
    print("\nVerdicts by source:")
    print(pd.crosstab(df["source"], df["verdict"]).to_string())
    print("\nModels:")
    print(df["model"].value_counts().to_string())


if __name__ == "__main__":
    main()
