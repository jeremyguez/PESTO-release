"""Draw 50 genes at random that already have PEPPER, GO and LOEUF.

The 400-gene universe is 75% an autism list. These fifty are a dilution
control: protein-coding, annotated on the three figure-3 tracks, and not
already scored. Six phenotypes, 300 pairs.

  python3 scripts/92_fig3_random50.py --write-pairs
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig3_random50/pairs.tsv --workers 20 \
      --runs-dir results/fig3_random50/runs/flow \
      --out results/fig3_random50/pesto.json
  python3 scripts/92_fig3_random50.py --assemble

--assemble concatenates the draw onto results/fig3_universe400_scores.tsv and
writes results/fig3_universe450_scores.tsv, which is what every downstream
figure should read from here on.

Usage: python3 scripts/92_fig3_random50.py --write-pairs | --assemble
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
# PEPPER gene features (Guez et al. 2026) are not redistributed here; place
# gene_features_for_s_het.tsv.gz and obs_exp_for_loeuf_missense.tsv in this
# directory to draw the random genes again.
PEPPER = ROOT / "data" / "pepper"
LOEUF = ROOT / "data" / "loeuf_scores.tsv"
GMT = (
    ROOT / "data" / "go" / "GO_Biological_Process_2026.gmt",
    ROOT / "data" / "go" / "GO_Molecular_Function_2026.gmt",
    ROOT / "data" / "go" / "GO_Cellular_Component_2026.gmt",
    ROOT / "data" / "syngo" / "SynGO_2024.gmt",
)
SCORED = ROOT / "results" / "fig3_universe400_scores.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe450_scores.tsv"
OUT = ROOT / "results" / "fig3_random50"
GENES = OUT / "genes.tsv"
PAIRS = OUT / "pairs.tsv"
PESTO = OUT / "pesto.tsv"

PHENOS = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
]
N = 50
SEED = 20260817


def annotated():
    mapping = pd.read_csv(
        PEPPER / "obs_exp_for_loeuf_missense.tsv", sep="\t",
        usecols=["gene_symbol", "ensg"])
    mapping = mapping[mapping.gene_symbol.notna() & mapping.ensg.notna()]
    mapping["gene"] = mapping.gene_symbol.str.upper()
    feat = pd.read_csv(
        PEPPER / "gene_features_for_s_het.tsv.gz",
        sep="\t", compression="gzip")
    pep = set(mapping.loc[mapping.ensg.isin(set(feat.ensg.dropna())), "gene"])

    loeuf = pd.read_csv(LOEUF, sep="\t")
    loeuf["gene"] = loeuf.gene_symbol.astype(str).str.upper()
    loeuf["LOEUF"] = pd.to_numeric(loeuf.LOEUF, errors="coerce")
    loe = set(loeuf.loc[loeuf.LOEUF.notna(), "gene"])

    go = set()
    for path in GMT:
        with open(path) as fh:
            for line in fh:
                go.update(g.upper() for g in line.rstrip("\n").split("\t")[2:] if g)
    return pep & loe & go


def already_scored():
    raw = pd.read_csv(SCORED, sep="\t")
    return set(raw.gene.str.upper())


def write_pairs():
    pool = sorted(annotated() - already_scored())
    rng = np.random.default_rng(SEED)
    pick = sorted(rng.choice(pool, N, replace=False).tolist())
    OUT.mkdir(parents=True, exist_ok=True)
    with GENES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "set"], delimiter="\t")
        w.writeheader()
        w.writerows({"gene": g, "set": "random50_ann"} for g in pick)
    rows = [{"gene": g, "phenotype": p, "set": "random50_ann"}
            for g in pick for p in PHENOS]
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
    new = new[["gene", "phenotype", "lit_mean"]].copy()
    new["source"] = "random50_ann"
    new["set"] = "random50_ann"

    clash = set(old.gene) & set(new.gene)
    if clash:
        raise SystemExit(f"draw overlaps the scored universe: {sorted(clash)}")

    out = pd.concat([old, new], ignore_index=True)
    out = out.sort_values(["gene", "phenotype"])
    out.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(out)} rows  "
          f"{out.gene.nunique()} genes")
    for p in PHENOS:
        sub = out[out.phenotype == p]
        reach = sub.lit_mean.isin(("Existing", "Established")).sum()
        print(f"  {p:<28} {sub.gene.nunique():4} scored  {reach:4} reaching")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-pairs", action="store_true")
    ap.add_argument("--assemble", action="store_true")
    args = ap.parse_args()
    if args.write_pairs:
        write_pairs()
    if args.assemble:
        assemble()
    if not args.write_pairs and not args.assemble:
        ap.error("pass --write-pairs and/or --assemble")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
