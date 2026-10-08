#!/usr/bin/env python3
"""Supplementary Table 1: every RVAS association of Figure 1 and its verdicts.

One row per unique gene-phenotype pair (950), read from the table behind
figure1_v7. A pair significant in both studies is listed once, with study
"Both" and the All of Us-GeneBass P value, as the table already keeps it.
Answers the request that the Novel and Hypothesized lists be released.

Writes results/supplementary_table1.tsv and, when openpyxl is installed,
results/supplementary_table1.xlsx with a column legend on a second sheet.

  python3 scripts/145_supplementary_table1.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "results" / "all_runs_auto_v7.tsv"
BRAVA_RAW = ROOT / "data" / "raw" / "Duncan_results.tsv"
OUT = ROOT / "results" / "supplementary_table1"

VERDICTS = ["Novel", "Hypothesized", "Existing", "Established"]
RANK = {v: i for i, v in enumerate(VERDICTS)}
PHENOTYPE_ALIASES = {"bmi": "body mass index"}

LEGEND = [
    ("gene", "HGNC gene symbol as queried."),
    ("phenotype", "Phenotype as named in the source study."),
    ("study", "All of Us-GeneBass meta-analysis, BRaVa, or Both."),
    ("p_value", "All of Us-GeneBass burden meta-analysis P; minimum P reported by "
                "BRaVa for BRaVa-only pairs."),
    ("literature_verdict", "Literature branch verdict: rounded confidence-weighted mean."),
    ("w_novel", "Literature branch confidence points on Novel (of 100)."),
    ("w_hypothesized", "Confidence points on Hypothesized."),
    ("w_existing", "Confidence points on Existing."),
    ("w_established", "Confidence points on Established."),
    ("evidence_score", "Confidence-weighted mean on the 1-4 scale."),
    ("n_pubmed_records", "PubMed records retrieved."),
    ("n_records_read", "Records passed to the final assessment."),
    ("cited_pmids", "PMIDs cited in the justification."),
    ("literature_justification", "Justification returned by the literature branch."),
    ("open_targets_verdict", "Open Targets branch verdict."),
    ("open_targets_trait", "Open Targets trait the verdict rests on, if any."),
    ("open_targets_score", "Open Targets association score of that trait."),
    ("max_verdict", "Stronger of the two branch verdicts."),
]


def norm_ph(x) -> str:
    s = re.sub(r"[^a-z0-9]+", " ", str(x).lower()).strip()
    return PHENOTYPE_ALIASES.get(s, s)


def main():
    d = pd.read_csv(RUNS, sep="\t")
    d = d[d.source.isin(["AoU", "BRAVA"]) & d.verdict.isin(VERDICTS)
          & d.ot_verdict.isin(VERDICTS)].copy()

    raw = pd.read_csv(BRAVA_RAW, sep="\t")
    brava = (set(zip(raw.external_gene_name.str.upper(), raw.phenotype.map(norm_ph)))
             | set(zip(raw.external_gene_name.str.upper(), raw.phenotype_full.map(norm_ph))))
    key = list(zip(d.gene.str.upper(), d.phenotype.map(norm_ph)))
    shared = np.array([k in brava for k in key]) & (d.source == "AoU").to_numpy()
    study = np.where(shared, "Both",
                     np.where(d.source == "AoU", "All of Us-GeneBass", "BRaVa"))

    w = d[["p_novel", "p_hypothesized", "p_existing", "p_established"]].astype(float)
    score = (w.to_numpy() * np.array([1, 2, 3, 4])).sum(axis=1) / w.sum(axis=1).to_numpy()
    mx = np.maximum(d.verdict.map(RANK), d.ot_verdict.map(RANK))

    t = pd.DataFrame({
        "gene": d.gene,
        "phenotype": d.phenotype,
        "study": study,
        "p_value": pd.to_numeric(d.burden_pvalue, errors="coerce"),
        "literature_verdict": d.verdict,
        "w_novel": w.p_novel.astype(int),
        "w_hypothesized": w.p_hypothesized.astype(int),
        "w_existing": w.p_existing.astype(int),
        "w_established": w.p_established.astype(int),
        "evidence_score": score.round(2),
        "n_pubmed_records": d.n_pubmed,
        "n_records_read": d.n_read,
        "cited_pmids": d.cited_pmids.fillna(""),
        "literature_justification": d.justification.fillna(""),
        "open_targets_verdict": d.ot_verdict,
        "open_targets_trait": d.ot_matched_traits.fillna(""),
        "open_targets_score": pd.to_numeric(d.ot_max_score, errors="coerce").round(3),
        "max_verdict": [VERDICTS[i] for i in mx],
    })
    t["_rank"] = t.max_verdict.map(RANK)
    t = t.sort_values(["_rank", "p_value", "gene"]).drop(columns="_rank")

    t.to_csv(OUT.with_suffix(".tsv"), sep="\t", index=False)
    print(f"wrote {OUT.with_suffix('.tsv').relative_to(ROOT)}: {len(t)} pairs")
    print(t.study.value_counts().to_dict())
    print(t.max_verdict.value_counts().reindex(VERDICTS).to_dict())
    try:
        with pd.ExcelWriter(OUT.with_suffix(".xlsx")) as xw:
            t.to_excel(xw, sheet_name="Supplementary Table 1", index=False)
            pd.DataFrame(LEGEND, columns=["column", "description"]).to_excel(
                xw, sheet_name="Columns", index=False)
        print(f"wrote {OUT.with_suffix('.xlsx').relative_to(ROOT)}")
    except ImportError:
        print("openpyxl not installed: TSV only")


if __name__ == "__main__":
    main()
