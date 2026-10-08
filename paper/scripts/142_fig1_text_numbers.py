#!/usr/bin/env python3
"""Every number the main text and legends quote about Figure 1 and ED Fig. 1.

Answers: how many association-study records and unique gene-phenotype pairs
the RVAS application covers, and what the literature, Open Targets and Max
verdicts, the negative-control contrast, the branch concordance and the
burden-significance trend come to on them. Reads the same two tables as
figure1_v7 (results/all_runs_auto_v7.tsv and the v7 absent-50 controls), so
the text and the figure cannot drift apart.

Shared pairs are counted once, on the AoU side, as in the table; the BRaVa
total adds them back from data/raw/Duncan_results.tsv.

  python3 scripts/142_fig1_text_numbers.py > results/fig1_text_numbers.txt
"""
from __future__ import annotations

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "results", "all_runs_auto_v7.tsv")
ABSENT = os.path.join(ROOT, "results", "bench_fig1c_absent50_hgnc_pesto_v7.tsv")
BRAVA_RAW = os.path.join(ROOT, "data", "raw", "Duncan_results.tsv")

VERDICTS = ["Novel", "Hypothesized", "Existing", "Established"]
RANK = {v: i for i, v in enumerate(VERDICTS)}
PHENOTYPE_ALIASES = {"bmi": "body mass index"}
MIN_PAIRS = 8

# Mirrors canon_pheno in scripts/25_figure_ed1.R.
ED1_CANON = {
    "height": "Height",
    "total cholesterol": "Total cholesterol",
    "waist to hip ratio adjusted for bmi": "WHRadjBMI",
    "hip-circumference-mean": "Hip circumference",
    "mean corpuscular volume": "Mean corpuscular volume",
    "red cell distribution width": "Red cell distribution width",
    "mean platelet volume": "Mean platelet volume",
    "bmi": "Body mass index",
}


def norm_ph(x) -> str:
    s = re.sub(r"[^a-z0-9]+", " ", str(x).lower()).strip()
    return PHENOTYPE_ALIASES.get(s, s)


def pct(k, n) -> str:
    return f"{k} ({100 * k / n:.0f}%)"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", default=RUNS)
    ap.add_argument("--absent", default=ABSENT)
    args = ap.parse_args()

    d = pd.read_csv(args.runs, sep="\t")
    d = d[d.source.isin(["AoU", "BRAVA"]) & d.verdict.isin(VERDICTS)
          & d.ot_verdict.isin(VERDICTS)].copy()
    d["key"] = list(zip(d.gene.str.upper(), d.phenotype.map(norm_ph)))
    d["lit"] = d.verdict.map(RANK)
    d["ot"] = d.ot_verdict.map(RANK)
    d["mx"] = np.maximum(d.lit, d.ot)
    n = len(d)

    raw = pd.read_csv(BRAVA_RAW, sep="\t")
    brava_keys = (set(zip(raw.external_gene_name.str.upper(), raw.phenotype.map(norm_ph)))
                  | set(zip(raw.external_gene_name.str.upper(), raw.phenotype_full.map(norm_ph))))
    aou = d[d.source == "AoU"]
    shared = int(sum(k in brava_keys for k in aou.key))
    n_aou = len(aou)
    n_brava_rows = int((d.source == "BRAVA").sum())
    n_brava = n_brava_rows + shared

    print("== Records and pairs")
    print(f"table rows (analysed)            {n}")
    print(f"unique gene-phenotype pairs      {len(set(d.key))}")
    print(f"AoU + GeneBass                   {n_aou}")
    print(f"BRaVa rows in table              {n_brava_rows}")
    print(f"AoU pairs also in BRaVa          {shared}")
    print(f"BRaVa evaluated (rows + shared)  {n_brava}  (raw BRaVa rows: {len(raw)})")
    print(f"association-study records        {n_aou + n_brava}")

    print("\n== Verdicts on the unique pairs")
    print(f"literature Novel                 {pct(int((d.lit == 0).sum()), n)}")
    print(f"Open Targets Novel               {pct(int((d.ot == 0).sum()), n)}")
    max_novel, max_hyp = int((d.mx == 0).sum()), int((d.mx == 1).sum())
    print(f"Max Novel                        {pct(max_novel, n)}")
    print(f"Max Hypothesized                 {pct(max_hyp, n)}")
    print(f"Max Novel or Hypothesized        {pct(max_novel + max_hyp, n)}")
    novel = d[d.lit == 0]
    n_pub = pd.to_numeric(novel.n_pubmed, errors="coerce")
    print(f"literature Novel with 0 records  {int((n_pub == 0).sum())} of {len(novel)}"
          f" ({100 * (n_pub == 0).sum() / len(novel):.0f}%), median records {n_pub.median():.0f}")

    c = pd.read_csv(args.absent, sep="\t")
    c = c[c.verdict.isin(VERDICTS) & c.open_targets_verdict.isin(VERDICTS)]
    c_mx = np.maximum(c.verdict.map(RANK), c.open_targets_verdict.map(RANK))
    c_nh = int((c_mx <= 1).sum())
    rv_nh = max_novel + max_hyp
    _, p_fisher = fisher_exact([[c_nh, len(c) - c_nh], [rv_nh, n - rv_nh]])
    print("\n== Negative controls")
    print(f"controls Max Novel or Hyp.       {c_nh} of {len(c)} ({100 * c_nh / len(c):.0f}%)")
    print(f"Fisher exact P vs RVAS           {p_fisher:.2g}")

    print("\n== Branch concordance")
    print(f"exact agreement                  {100 * (d.lit == d.ot).mean():.1f}%")
    print(f"within one level                 {100 * ((d.lit - d.ot).abs() <= 1).mean():.1f}%")

    p = pd.to_numeric(d.burden_pvalue, errors="coerce")
    ok = p.notna() & (p > 0)
    logp = -np.log10(p[ok])
    rho, p_rho = spearmanr(d.mx[ok], logp)
    print("\n== Burden significance by Max verdict")
    print(f"pairs with a P value             {int(ok.sum())}")
    print(f"Spearman rho                     {rho:.2f} (P = {p_rho:.2g})")
    for i, v in enumerate(VERDICTS):
        sel = d.mx[ok] == i
        print(f"mean -log10 P {v:<13}      {logp[sel].mean():.1f}  (n = {int(sel.sum())})")

    ph = d.phenotype.map(lambda x: ED1_CANON.get(str(x).strip().lower(), str(x).strip()))
    counts = ph.value_counts()
    keep = counts[counts >= MIN_PAIRS]
    print("\n== ED Fig. 1b,c")
    print(f"phenotypes with >= {MIN_PAIRS} pairs        {len(keep)} covering {int(keep.sum())} pairs")


if __name__ == "__main__":
    main()
