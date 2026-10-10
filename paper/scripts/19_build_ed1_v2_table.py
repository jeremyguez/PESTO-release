#!/usr/bin/env python3
"""Build a second Extended Data Table 1, grouped by evidence class.

Eight pairs chosen to show the kinds of support a Max Novel or Hypothesized
call can rest on: no report at all, a neighbouring phenotype, a cell mechanism,
an animal model, an indirect human finding, or an Open Targets tie that is not
the queried trait. The original ten-row table is left untouched.

Numbers come from results/all_runs_auto_v7.tsv so the verdicts match Figure 1
v7. A pair whose literature or Open Targets verdict has risen past
Hypothesized aborts the build.

Usage: python3 scripts/19_build_ed1_v2_table.py
"""

import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "results", "all_runs_auto_v7.tsv")
ALIASES = os.path.join(ROOT, "data", "phenotype_aliases.tsv")
OUT = os.path.join(ROOT, "results", "ed1_v2_pairs.tsv")

ALLOWED = {"Novel", "Hypothesized"}

# Order is the argument, not burden P.
CURATION = [
    {
        "gene": "THNSL1",
        "phenotype": "Stroke",
        "category": "Only Novel (insufficient)",
        "evidence": "NA",
        "interest": (
            "THNSL1 has not been directly implicated in stroke; it was only "
            "among >1,900 genes altered in mouse splenic B cells after "
            "protective hypoxic preconditioning, before stroke induction, "
            "with no functional follow-up."
        ),
        "refs": "Monson NL et al. 2014",
        "pmids": "24485041",
    },
    {
        "gene": "NAA15",
        "phenotype": "Hypertension",
        "category": "Only Novel (insufficient)",
        "evidence": "NA",
        "interest": (
            "NAA15 haploinsufficiency causes congenital heart disease, but "
            "hypertension was reported in only 1/19 evaluated individuals "
            "ascertained for truncating variants (Cheng H et al. 2018)."
        ),
        "refs": "Cheng H et al. 2018; Ward T et al. 2021",
        "pmids": "29656860;33557580",
    },
    {
        "gene": "OSBP",
        "phenotype": "HDL cholesterol",
        "category": "Hypothesized, mechanistic",
        "evidence": "mechanistic",
        "interest": (
            "OSBP knockdown increases cholesterol efflux in CHO cells "
            "and raises ABCA1 in J774 macrophages, providing a "
            "mechanistic link to HDL biology without direct human genetic "
            "evidence."
        ),
        "refs": "Bowden K, Ridgway ND. 2008",
        "pmids": "18450749",
    },
    {
        "gene": "PODN",
        "phenotype": "Coronary artery disease",
        "category": "Hypothesized, animal",
        "evidence": "animal / functional",
        "interest": (
            "Podn-null mice show excessive neointimal lesion formation "
            "after arterial injury; overexpression of human PODN "
            "inhibits human vascular smooth-muscle-cell migration and "
            "proliferation, and PODN is expressed in human atheroma."
        ),
        "refs": "Hutter R et al. 2013",
        "pmids": "24043300",
    },
    {
        "gene": "SIN3A",
        "phenotype": "Stroke",
        "category": "Hypothesized on both",
        "evidence": "mechanistic / animal",
        "interest": (
            "Cerebral ischemia recruits mSin3A to REST-corepressor "
            "complexes, providing mechanistic support without direct "
            "human genetic evidence for stroke."
        ),
        "refs": "Noh KM et al. 2012",
        "pmids": "22371606",
    },
    {
        "gene": "IRS2",
        "phenotype": "Acute kidney failure",
        "category": "Hypothesized on both",
        "evidence": "indirect human genetic + animal",
        "interest": (
            "IRS2 truncating variants increase chronic kidney disease "
            "risk independently of diabetes, and ghrelin-mediated "
            "protection from ischemic acute kidney injury is lost in "
            "Irs2-knockout mice."
        ),
        "refs": "Zhao Y et al. 2025; Takeda R et al. 2006",
        "pmids": "41073786;16306169",
    },
    {
        "gene": "CPXM2",
        "phenotype": "Rheumatoid arthritis",
        "category": "Novel lit, Hypoth. OT",
        "evidence": "transcriptomic / related disease",
        "interest": (
            "No retrieved literature directly links CPXM2 to rheumatoid "
            "arthritis; Open Targets provides low-confidence "
            "transcriptomic evidence of differential CPXM2 expression in "
            "rheumatoid arthritis and additional indirect support through "
            "ulcerative colitis, a related immune-mediated inflammatory "
            "disease."
        ),
        "refs": "-",
        "pmids": "",
    },
    {
        "gene": "SETD1A",
        "phenotype": "Asthma",
        "category": "Novel lit, Hypoth. OT",
        "evidence": "lung function",
        "interest": (
            "No retrieved literature directly links SETD1A to asthma. "
            "Open Targets provides indirect support through forced "
            "expiratory volume, with a low evidence score (0.053)."
        ),
        "refs": "-",
        "pmids": "",
    },
]

FIELDS = [
    "gene", "phenotype", "cohort", "burden_pvalue", "p_fmt",
    "lit_verdict", "ot_verdict", "category", "evidence",
    "n_pubmed", "interest", "refs", "pmids",
]


with open(ALIASES, encoding="utf-8") as _fh:
    PHENOTYPE_ALIASES = {r["alias"]: r["canonical"]
                         for r in csv.DictReader(_fh, delimiter="\t")}


def pair_key(gene, phenotype):
    # A pair significant in both studies keeps its AoU row, which may name
    # the phenotype differently from the curation (HDL cholesterol).
    ph = re.sub(r"[^a-z0-9]+", " ", phenotype.lower()).strip()
    return gene.upper(), PHENOTYPE_ALIASES.get(ph, ph)


def load(path):
    if not os.path.exists(path):
        sys.exit(f"missing input: {path}")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {pair_key(r["gene"], r["phenotype"]): r for r in rows}


def main():
    runs = load(RUNS)
    out, problems = [], []
    for c in CURATION:
        key = pair_key(c["gene"], c["phenotype"])
        r = runs.get(key)
        if r is None:
            problems.append(f"{c['gene']} / {c['phenotype']}: absent from {os.path.basename(RUNS)}")
            continue
        lit, otv = r["verdict"], r["ot_verdict"]
        if lit not in ALLOWED or otv not in ALLOWED:
            problems.append(
                f"{c['gene']} / {c['phenotype']}: verdicts are {lit} (literature) "
                f"and {otv} (Open Targets); one has passed Hypothesized"
            )
            continue
        p = float(r["burden_pvalue"])
        out.append({
            "gene": c["gene"],
            "phenotype": c["phenotype"],
            "cohort": r["source"],
            "burden_pvalue": r["burden_pvalue"],
            "p_fmt": f"{p:.1e}",
            "lit_verdict": lit,
            "ot_verdict": otv,
            "category": c["category"],
            "evidence": c["evidence"],
            "n_pubmed": r.get("n_pubmed", ""),
            "interest": c["interest"],
            "refs": c["refs"],
            "pmids": c["pmids"],
        })

    if problems:
        sys.exit("curation no longer matches the runs:\n  " + "\n  ".join(problems))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"Wrote {OUT} ({len(out)} pairs)")
    for r in out:
        print(f"  {r['category'][:28]:28} {r['gene']:8} {r['phenotype'][:32]:32} "
              f"{r['cohort']:6} {r['p_fmt']:>8}  lit={r['lit_verdict']:12} "
              f"ot={r['ot_verdict']:12} ev={r['evidence']}")


if __name__ == "__main__":
    main()
