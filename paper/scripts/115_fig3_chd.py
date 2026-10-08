"""Congenital heart disease on the 509-gene Alzheimer universe.

Adds the 60 PCGC-significant dominant CHD genes (Sierant 2025, Dataset S7).
Scores congenital heart disease on every gene already in the 509 table, and
scores the genes that are new to that table on the seven existing phenotypes.
Arm `current`, no Open Targets. Does not overwrite the 509-gene files or
figures/figure3.*

  python3 scripts/115_fig3_chd.py --write-pairs
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig3_chd/pairs.tsv --workers 20 \
      --runs-dir results/fig3_chd/runs/flow \
      --out results/fig3_chd/pesto.json
  python3 scripts/115_fig3_chd.py --assemble

Usage: python3 scripts/115_fig3_chd.py --write-pairs | --assemble
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "matrix5", ROOT / "scripts" / "86_go_matrix_5pheno.py")
g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(g)

OLD = ROOT / "results" / "fig3_universe509_ad_scores.tsv"
GENES = ROOT / "data" / "chd" / "pcgc_sig60.tsv"
OUT = ROOT / "results" / "fig3_chd"
PAIRS = OUT / "pairs.tsv"
NEW_GENES = OUT / "genes_new.tsv"
PESTO = OUT / "pesto.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe551_chd_scores.tsv"

CHD = "congenital heart disease"
OTHER = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
    "Alzheimer's disease",
]
SET_OLD = "universe509"
SET_NEW = "pcgc_sig60"


def pcgc60():
    return sorted(pd.read_csv(GENES, sep="\t").gene.str.upper().unique())


def write_pairs():
    old = pd.read_csv(OLD, sep="\t")
    old["gene"] = old.gene.str.upper()
    have = set(old.gene)
    sig = pcgc60()
    new = [g for g in sig if g not in have]
    if len(sig) != 60:
        raise SystemExit(f"expected 60 PCGC genes, got {len(sig)}")
    universe = sorted(have | set(sig))
    rows = [{"gene": g, "phenotype": CHD, "set": SET_OLD if g in have else SET_NEW}
            for g in universe]
    rows += [{"gene": g, "phenotype": p, "set": SET_NEW}
             for g in new for p in OTHER]
    OUT.mkdir(parents=True, exist_ok=True)
    with NEW_GENES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "set"], delimiter="\t")
        w.writeheader()
        w.writerows({"gene": g, "set": SET_NEW} for g in new)
    with PAIRS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"universe {len(have)} + {len(new)} new = {len(universe)}")
    print(f"already in 509: {len(sig) - len(new)}")
    print(f"wrote {NEW_GENES.relative_to(ROOT)}")
    print(f"wrote {PAIRS.relative_to(ROOT)}  {len(rows)} pairs  "
          f"{len(universe)} x {CHD!r} + {len(new)} x {len(OTHER)} other")


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
    keep["source"] = "fig3_chd"
    keep["set"] = [set_of.get((r.gene, r.phenotype), "")
                   for r in keep.itertuples()]
    out = pd.concat([old, keep], ignore_index=True)
    out = out.drop_duplicates(["gene", "phenotype"], keep="last")
    out = out.sort_values(["gene", "phenotype"])
    out.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(out)} rows  "
          f"{out.gene.nunique()} genes")
    for p in OTHER + [CHD]:
        sub = out[out.phenotype == p]
        reach = sub.lit_mean.isin(("Existing", "Established")).sum()
        print(f"  {p:<28} {sub.gene.nunique():4} scored  {reach:4} reaching")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-pairs", action="store_true")
    ap.add_argument("--assemble", action="store_true")
    args = ap.parse_args()
    if args.write_pairs:
        write_pairs()
    if args.assemble:
        assemble()
    if not (args.write_pairs or args.assemble):
        ap.error("pass --write-pairs and/or --assemble")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
