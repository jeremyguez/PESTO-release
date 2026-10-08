"""Alzheimer's disease on the figure-3 universe, plus nine GenCC genes.

Scores the 450 genes of fig3_universe450_scores.tsv and nine GenCC Alzheimer
genes (A2M and PRNP dropped) for Alzheimer's disease, and scores those nine
on the six figure-3 phenotypes they lack. Arm `current`, no Open Targets.
Writes a parallel scores table and GO inputs; does not touch the 450-gene
files or figures/figure3.*

  python3 scripts/104_fig3_alzheimer.py --write-pairs
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig3_alzheimer/pairs.tsv --workers 20 \
      --runs-dir results/fig3_alzheimer/runs/flow \
      --out results/fig3_alzheimer/pesto.json
  python3 scripts/104_fig3_alzheimer.py --assemble --go
  Rscript scripts/105_figure3_ad.R
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "matrix5", ROOT / "scripts" / "86_go_matrix_5pheno.py")
g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(g)

OLD = ROOT / "results" / "fig3_universe450_scores.tsv"
OUT = ROOT / "results" / "fig3_alzheimer"
PAIRS = OUT / "pairs.tsv"
PESTO = OUT / "pesto.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe459_ad_scores.tsv"
SIG = OUT / "go_signatures.tsv"
CORR = OUT / "pheno_go_correlation.tsv"
HEAT = OUT / "pheno_heatmap.tsv"

AD = "Alzheimer's disease"
OTHER = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
]
GENCC9 = [
    "ABCA7", "APOE", "APP", "ECE2", "GRIN2C",
    "PSEN1", "PSEN2", "SORL1", "TREM2",
]

# Signature pick matches scripts/93_go_signatures450.py (DD omitted there).
SIG_PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
    ("epilepsy", "epilepsy"),
    ("bipolar disorder", "bipolar"),
    ("schizophrenia", "SCZ"),
    ("type 2 diabetes", "T2D"),
    (AD, "AD"),
]
# Clustering matches scripts/96_pheno_dendrogram.py plus AD.
CLUST_PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
    ("developmental disorder", "DD"),
    ("epilepsy", "epilepsy"),
    ("bipolar disorder", "bipolar"),
    ("schizophrenia", "SCZ"),
    ("type 2 diabetes", "T2D"),
    (AD, "AD"),
]
FDR = 0.05
OVERLAP = 0.5
TOP_PER_PHENO = 6
MIN_K_PICK = 5
MIN_K_CLUST = 10


def write_pairs():
    genes450 = sorted(set(pd.read_csv(OLD, sep="\t").gene.str.upper()))
    new9 = [g for g in GENCC9 if g not in set(genes450)]
    if len(new9) != 9:
        raise SystemExit(f"expected 9 new genes, got {new9}")
    rows = [{"gene": g, "phenotype": AD, "set": "universe450"}
            for g in genes450]
    rows += [{"gene": g, "phenotype": AD, "set": "gencc9_ad"}
             for g in new9]
    rows += [{"gene": g, "phenotype": p, "set": "gencc9_ad"}
             for g in new9 for p in OTHER]
    OUT.mkdir(parents=True, exist_ok=True)
    with PAIRS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {PAIRS.relative_to(ROOT)}  {len(rows)} pairs  "
          f"{len(genes450)} x {AD!r} + {len(new9)} x "
          f"{AD!r} + {len(new9)} x {len(OTHER)} other")


def assemble():
    old = pd.read_csv(OLD, sep="\t")
    old["gene"] = old.gene.str.upper()
    new = pd.read_csv(PESTO, sep="\t").rename(columns={"verdict": "lit_mean"})
    new["gene"] = new.gene.str.upper()
    pairs = pd.read_csv(PAIRS, sep="\t")
    pairs["gene"] = pairs.gene.str.upper()
    set_of = {(r.gene, r.phenotype): r.set for r in pairs.itertuples()}
    keep = new[new.lit_mean.isin(g.W)].copy()
    missing = new[~new.lit_mean.isin(g.W)]
    if len(missing):
        print(f"dropping {len(missing)} rows without a lit_mean:")
        print(missing[["gene", "phenotype", "lit_mean", "error"]]
              .to_string(index=False))
    keep = keep[["gene", "phenotype", "lit_mean"]].copy()
    keep["source"] = "fig3_alzheimer"
    keep["set"] = [set_of.get((r.gene, r.phenotype), "")
                   for r in keep.itertuples()]
    out = pd.concat([old, keep], ignore_index=True)
    out = out.drop_duplicates(["gene", "phenotype"], keep="last")
    out = out.sort_values(["gene", "phenotype"])
    out.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(out)} rows  "
          f"{out.gene.nunique()} genes")
    for p in OTHER + [AD]:
        sub = out[out.phenotype == p]
        reach = sub.lit_mean.isin(("Existing", "Established")).sum()
        print(f"  {p:<28} {sub.gene.nunique():4} scored  {reach:4} reaching")


def pick(by_item):
    cand = [r for r in by_item.values()
            if r["q"] < FDR and r["k"] >= MIN_K_PICK]
    cand.sort(key=lambda r: (r["q"], -r["lor"]))
    kept = []
    for r in cand:
        if any(len(r["hit"] & s["hit"]) / min(len(r["hit"]), len(s["hit"]))
               >= OVERLAP for s in kept):
            continue
        kept.append(r)
        if len(kept) == TOP_PER_PHENO:
            break
    return kept


def go():
    raw = pd.read_csv(UNIVERSE, sep="\t")
    raw = raw[raw.lit_mean.isin(g.W)]
    raw["gene"] = raw.gene.str.upper()
    n_genes = raw.gene.nunique()

    sig_tables, chosen = {}, {}
    for pheno, short in SIG_PHENOS:
        sub = raw[raw.phenotype == pheno]
        score = {gene: g.W[v] for gene, v in zip(sub.gene, sub.lit_mean)}
        sig_tables[short] = enrich = g.enrich(score)
        chosen[short] = pick(enrich)
        n_fg = next(iter(enrich.values()))["n_fg"]
        n_sig = sum(r["q"] < FDR and r["k"] >= MIN_K_PICK
                    for r in enrich.values())
        head = chosen[short][0]["item"] if chosen[short] else "-- none --"
        print(f"{short:8} reaching {n_fg:3}/{n_genes}  FDR<{FDR}: {n_sig:3}"
              f"  picked {len(chosen[short])}  top: {head[:44]}")

    selected = {r["item"] for picks in chosen.values() for r in picks}
    order = {short: i for i, (_, short) in enumerate(SIG_PHENOS)}
    block = {}
    for item in selected:
        owners = [(sig_tables[s].get(item), s) for _, s in SIG_PHENOS]
        owners = [(r, s) for r, s in owners if r is not None and r["q"] < FDR]
        best = max(owners, key=lambda rs: rs[0]["lor"])
        block[item] = (best[1], best[0]["q"])

    rows = []
    for item in sorted(selected, key=lambda i: (order[block[i][0]], block[i][1])):
        source, _ = block[item]
        for _, short in SIG_PHENOS:
            r = sig_tables[short].get(item)
            if r is None:
                continue
            rows.append({
                "item": item, "source": source, "phenotype": short,
                "k": r["k"], "K": r["K"], "n_fg": r["n_fg"],
                "lor": r["lor"], "q": r["q"], "p": r["p"],
                "frac": r["k"] / r["K"],
            })
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(SIG, sep="\t", index=False)
    print(f"wrote {SIG.relative_to(ROOT)}  terms={len(selected)}")

    tables = {}
    for pheno, short in CLUST_PHENOS:
        sub = raw[raw.phenotype == pheno]
        score = {gene: g.W[v] for gene, v in zip(sub.gene, sub.lit_mean)}
        tables[short] = enrich = g.enrich(score)
        n_fg = next(iter(enrich.values()))["n_fg"]
        n_sig = sum(r["q"] < FDR for r in enrich.values())
        print(f"{short:8} cluster reaching {n_fg:3}/{n_genes}  "
              f"FDR<{FDR}: {n_sig:3}")

    shared = set.intersection(*(set(t) for t in tables.values()))
    shared = {i for i in shared
              if next(iter(tables.values()))[i]["K"] >= MIN_K_CLUST}
    items = sorted(shared)
    shorts = [s for _, s in CLUST_PHENOS]
    mat = np.array([[tables[s][i]["lor"] for i in items] for s in shorts])
    rho = spearmanr(mat, axis=1).statistic
    corr = pd.DataFrame(rho, index=shorts, columns=shorts)
    corr.to_csv(CORR, sep="\t")
    print(f"\nclustering on {len(items)} terms with K >= {MIN_K_CLUST}")
    print(corr.round(3).to_string())

    dist = 1 - corr.values
    np.fill_diagonal(dist, 0.0)
    link = linkage(squareform(dist, checks=False), method="average")
    coph = cophenet(link, squareform(dist, checks=False))[0]
    leaves = dendrogram(link, no_plot=True, labels=shorts)["ivl"]
    print(f"\naverage-linkage order: {' '.join(leaves)}")
    print(f"cophenetic correlation: {coph:.3f}")

    shown = list(dict.fromkeys(pd.read_csv(SIG, sep="\t").item))
    heat = []
    for item in shown:
        for short in shorts:
            r = tables[short].get(item)
            if r is None:
                continue
            heat.append({
                "item": item, "phenotype": short, "lor": r["lor"],
                "q": r["q"], "k": r["k"], "K": r["K"], "n_fg": r["n_fg"],
            })
    pd.DataFrame(heat).to_csv(HEAT, sep="\t", index=False)
    print(f"wrote {CORR.relative_to(ROOT)}\nwrote {HEAT.relative_to(ROOT)}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-pairs", action="store_true")
    ap.add_argument("--assemble", action="store_true")
    ap.add_argument("--go", action="store_true")
    args = ap.parse_args()
    if args.write_pairs:
        write_pairs()
    if args.assemble:
        assemble()
    if args.go:
        go()
    if not (args.write_pairs or args.assemble or args.go):
        ap.error("pass --write-pairs, --assemble and/or --go")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
