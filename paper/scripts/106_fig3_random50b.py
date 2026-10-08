"""Second annotated random-50 draw, on the 459-gene Alzheimer universe.

Same pool as scripts/92_fig3_random50.py (PEPPER ∩ LOEUF ∩ GO), excluding
genes already in fig3_universe459_ad_scores.tsv. Seven phenotypes (the six
figure-3 traits plus Alzheimer's disease), 350 pairs. Does not touch the
450-gene files, the 459-gene table, or figures/figure3.*

  python3 scripts/106_fig3_random50b.py --write-pairs
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig3_random50b/pairs.tsv --workers 20 \
      --runs-dir results/fig3_random50b/runs/flow \
      --out results/fig3_random50b/pesto.json
  python3 scripts/106_fig3_random50b.py --assemble --go
  Rscript scripts/107_figure3_ad_r50.R

Usage: python3 scripts/106_fig3_random50b.py --write-pairs | --assemble | --go
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
SPEC92 = importlib.util.spec_from_file_location(
    "r92", ROOT / "scripts" / "92_fig3_random50.py")
r92 = importlib.util.module_from_spec(SPEC92)
SPEC92.loader.exec_module(r92)
SPEC = importlib.util.spec_from_file_location(
    "matrix5", ROOT / "scripts" / "86_go_matrix_5pheno.py")
g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(g)

SCORED = ROOT / "results" / "fig3_universe459_ad_scores.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe509_ad_scores.tsv"
OUT = ROOT / "results" / "fig3_random50b"
GENES = OUT / "genes.tsv"
PAIRS = OUT / "pairs.tsv"
PESTO = OUT / "pesto.tsv"
SIG = OUT / "go_signatures.tsv"
CORR = OUT / "pheno_go_correlation.tsv"
HEAT = OUT / "pheno_heatmap.tsv"

AD = "Alzheimer's disease"
PHENOS = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
    AD,
]
SET = "random50_ann2"
N = 50
SEED = 20260901

SIG_PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
    ("epilepsy", "epilepsy"),
    ("bipolar disorder", "bipolar"),
    ("schizophrenia", "SCZ"),
    ("type 2 diabetes", "T2D"),
    (AD, "AD"),
]
# DD dropped from this version: it tracks ASD and collapses the tree.
CLUST_PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
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


def already_scored():
    raw = pd.read_csv(SCORED, sep="\t")
    return set(raw.gene.str.upper())


def write_pairs():
    pool = sorted(r92.annotated() - already_scored())
    rng = np.random.default_rng(SEED)
    pick = sorted(rng.choice(pool, N, replace=False).tolist())
    OUT.mkdir(parents=True, exist_ok=True)
    with GENES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "set"], delimiter="\t")
        w.writeheader()
        w.writerows({"gene": gene, "set": SET} for gene in pick)
    rows = [{"gene": gene, "phenotype": p, "set": SET}
            for gene in pick for p in PHENOS]
    with PAIRS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"pool {len(pool)}  drew {len(pick)}  seed {SEED}")
    print(f"wrote {GENES.relative_to(ROOT)}")
    print(f"wrote {PAIRS.relative_to(ROOT)}  {len(rows)} pairs")
    print("genes: " + " ".join(pick))


def assemble():
    old = pd.read_csv(SCORED, sep="\t")
    old["gene"] = old.gene.str.upper()
    new = pd.read_csv(PESTO, sep="\t").rename(columns={"verdict": "lit_mean"})
    new["gene"] = new.gene.str.upper()
    keep = new[new.lit_mean.isin(g.W)].copy()
    missing = new[~new.lit_mean.isin(g.W)]
    if len(missing):
        print(f"dropping {len(missing)} rows without a lit_mean:")
        print(missing[["gene", "phenotype", "lit_mean", "error"]]
              .to_string(index=False))
    keep = keep[["gene", "phenotype", "lit_mean"]].copy()
    keep["source"] = SET
    keep["set"] = SET

    clash = set(old.gene) & set(keep.gene)
    if clash:
        raise SystemExit(f"draw overlaps the scored universe: {sorted(clash)}")

    out = pd.concat([old, keep], ignore_index=True)
    out = out.sort_values(["gene", "phenotype"])
    out.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(out)} rows  "
          f"{out.gene.nunique()} genes")
    for p in PHENOS:
        sub = out[out.phenotype == p]
        reach = sub.lit_mean.isin(("Existing", "Established")).sum()
        print(f"  {p:<28} {sub.gene.nunique():4} scored  {reach:4} reaching")


def pick(by_item, by="q"):
    cand = [r for r in by_item.values()
            if r["q"] < FDR and r["k"] >= MIN_K_PICK]
    if by == "or":
        cand.sort(key=lambda r: (-r["lor"], r["q"]))
    else:
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
        cation = [r for r in enrich.values() if "0098655" in r["item"]]
        if cation and short == "epilepsy":
            r = cation[0]
            print(f"         cation k={r['k']} K={r['K']} n_fg={r['n_fg']} "
                  f"p={r['p']:.4e} q={r['q']:.4f} lor={r['lor']:.2f}")

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
