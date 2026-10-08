"""Pairs still needed so figure 3 can use the 400-gene universe.

151 extras minus the 50 hardneg leaves 101 genes. They already have ASD
and SCZ; this writes the four missing phenotypes (404 pairs) and, after
scripts/run_arm.py has scored them, assembles results/fig3_universe400_scores.tsv.

  python3 scripts/82_fig3_four_phenos.py --write-pairs
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig3_four_phenos/pairs.tsv --workers 10 \
      --runs-dir results/fig3_four_phenos/runs/flow \
      --out results/fig3_four_phenos/pesto.json
  python3 scripts/82_fig3_four_phenos.py --assemble
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "fig3_four_phenos"
PAIRS = OUT / "pairs.tsv"
SCORED = OUT / "pesto.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe400_scores.tsv"

PHENOS = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
]
MISSING = [
    "developmental disorder",
    "epilepsy",
    "bipolar disorder",
    "type 2 diabetes",
]
BENCHES = (
    ("asc_ref299", ROOT / "benchmark" / "asc_ref299" / "pesto.tsv"),
    ("asd_extra101", ROOT / "benchmark" / "asd_extra101" / "pesto.tsv"),
    ("scz_extra101", ROOT / "benchmark" / "scz_extra101" / "pesto.tsv"),
    ("extra50_random", ROOT / "benchmark" / "extra50_random" / "pesto.tsv"),
)


def load(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def hardneg():
    return {r["gene"].upper() for r in load(BENCHES[1][1])
            if r.get("set") == "hardneg"}


def extras_keep():
    drop = hardneg()
    genes = {}
    for name, path in BENCHES[1:]:
        for r in load(path):
            g = r["gene"].upper()
            if g in drop:
                continue
            genes[g] = r.get("set") or ("random50" if name == "extra50_random"
                                        else "random")
    return genes


def write_pairs():
    genes = extras_keep()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [{"gene": g, "phenotype": p, "set": s}
            for g, s in sorted(genes.items()) for p in MISSING]
    with PAIRS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {PAIRS.relative_to(ROOT)}  {len(rows)} pairs  "
          f"{len(genes)} genes")


def assemble():
    drop = hardneg()
    rows = []
    seen = set()
    for name, path in BENCHES:
        for r in load(path):
            g = r["gene"].upper()
            if g in drop:
                continue
            p = r.get("phenotype")
            if p not in PHENOS:
                continue
            call = r.get("lit_mean") or r.get("verdict")
            if not call:
                continue
            key = (g, p)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"gene": g, "phenotype": p, "lit_mean": call,
                         "source": name, "set": r.get("set") or ""})
    if not SCORED.exists() and not (OUT / "pesto.tsv").exists():
        scored_path = ROOT / "results" / "arm_current_pairs.tsv"
    scored = []
    for path in (SCORED, OUT / "pesto.tsv",
                 ROOT / "results" / "fig3_four_phenos" / "pesto.tsv"):
        if path.exists():
            scored = load(path)
            break
    json_tsv = OUT / "pesto.tsv"
    # run_arm writes --out pesto.json and pesto.tsv next to it
    arm_tsv = OUT / "pesto.tsv"
    if not scored:
        cand = list(OUT.glob("*.tsv"))
        print("looking in", OUT, [p.name for p in cand])
    for r in scored:
        g = r["gene"].upper()
        p = r.get("phenotype")
        call = r.get("lit_mean") or r.get("verdict")
        if p not in PHENOS or not call or call == "":
            continue
        key = (g, p)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"gene": g, "phenotype": p, "lit_mean": call,
                     "source": "fig3_four_phenos", "set": r.get("set") or ""})
    UNIVERSE.parent.mkdir(parents=True, exist_ok=True)
    with UNIVERSE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "lit_mean", "source",
                                "set"], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["gene"], r["phenotype"])))
    genes = {r["gene"] for r in rows}
    by_p = {}
    for r in rows:
        by_p.setdefault(r["phenotype"], set()).add(r["gene"])
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(rows)} rows  "
          f"{len(genes)} genes")
    for p in PHENOS:
        print(f"  {p:<28} {len(by_p.get(p, ()))}")


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
