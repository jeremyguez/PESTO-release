"""Build one benchmark covering 299 genes on six phenotypes, reusing every run.

The 597 runs already paid for live in five separate benchmark folders, and the
flow cache is keyed by folder: a fresh benchmark would recompute all of them.
So this copies the existing run artefacts into the consolidated folder before
anything is launched, leaving only the genuinely missing pairs to compute.

Two normalisations differ and both must be respected, or the copy is wasted:

  flow/store.slug   lowercases, so gene and phenotype case do not matter
  the Open Targets cache filename keeps the case of the phenotype verbatim

That is why the phenotype strings below are copied character for character from
the source benchmarks -- "Autism Spectrum Disorder" is title case because the
autism benchmark wrote it that way, and lowercasing it would silently orphan
253 Open Targets files.

Reuse of the Open Targets files is only sound because those files carry no
version stamp and the matcher configuration has not moved; check that
config.ot_match_version() still reads tagged_shortlist:biolord:20:gate before
trusting them.
"""

import csv
import json
import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmark"
TARGET = BENCH / "asc_ref299"

SOURCES = ["asc_cluster2", "asc_cluster3", "asc_cluster5", "asc_cluster6", "asd"]

PHENOTYPES = [
    "Autism Spectrum Disorder",
    "developmental disorder",
    "epilepsy",
    "schizophrenia",
    "bipolar disorder",
    "type 2 diabetes",
]


def asc_genes():
    path = ROOT / "data" / "raw" / "ASC_GeneList_12082026.txt"
    with open(path, newline="") as fh:
        return {r["Gene"].strip().upper() for r in csv.DictReader(fh)
                if r["Gene"].strip()}


def reference_genes():
    sets = json.load(open(ROOT / "data" / "genelists" / "reference_sets.json"))
    schema = set(sets["schema_scz"])
    epi25 = set().union(*(set(v) for k, v in sets.items()
                          if k.startswith("epi25_")))
    return schema, epi25


def build_gene_list():
    asc = asc_genes()
    schema, epi25 = reference_genes()
    genes = sorted(asc | schema | epi25)
    print(f"ASC {len(asc)} | SCHEMA {len(schema)} (dont {len(schema & asc)} deja "
          f"dans ASC) | Epi25 {len(epi25)} (dont {len(epi25 & asc)} deja dans ASC)")
    print(f"union : {len(genes)} genes")
    return genes


def write_pairs(genes):
    TARGET.mkdir(parents=True, exist_ok=True)
    out = TARGET / "pairs.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "phenotype"])
        for g in genes:
            for p in PHENOTYPES:
                w.writerow([g, p])
    n = len(genes) * len(PHENOTYPES)
    print(f"ecrit {out} : {n} paires")
    return n


def copy_runs():
    flow_dir = TARGET / "runs" / "flow"
    flow_dir.mkdir(parents=True, exist_ok=True)
    n_flow = n_ot = skipped = 0
    for name in SOURCES:
        src = BENCH / name / "runs"
        for f in sorted((src / "flow").glob("*.json")):
            dst = flow_dir / f.name
            if dst.exists():
                skipped += 1
                continue
            shutil.copy2(f, dst)
            n_flow += 1
        for f in sorted(src.glob("open_targets_traits_*.txt")):
            dst = TARGET / "runs" / f.name
            if dst.exists():
                skipped += 1
                continue
            shutil.copy2(f, dst)
            n_ot += 1
    print(f"copie : {n_flow} runs flow, {n_ot} caches Open Targets"
          + (f", {skipped} deja presents" if skipped else ""))
    return n_flow, n_ot


def verify(total_pairs):
    """Count how many of the declared pairs the copied cache actually answers."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from pesto.flow import arms, store
    from pesto.flow.types import Query

    arm = arms.ARMS["current"]
    flow_dir = TARGET / "runs" / "flow"
    hits = misses = 0
    missing_ot = 0
    with open(TARGET / "pairs.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = Query(gene=row["gene"], phenotype=row["phenotype"])
            if Path(store.path(q, arm, runs_dir=str(flow_dir))).exists():
                hits += 1
                safe = row["phenotype"].replace(" ", "_").replace("/", "_")
                ot = TARGET / "runs" / f"open_targets_traits_{row['gene']}_{safe}.txt"
                if not ot.exists():
                    missing_ot += 1
            else:
                misses += 1
    print(f"\nverification sur {total_pairs} paires :")
    print(f"  deja repondues par le cache : {hits}")
    print(f"  a calculer                  : {misses}")
    if missing_ot:
        print(f"  ATTENTION {missing_ot} runs flow sans cache Open Targets associe")
    return hits, misses


def main():
    if (TARGET / "pairs.tsv").exists():
        print(f"{TARGET} existe deja ; rien n'est ecrase, les copies reprennent.")
    genes = build_gene_list()
    total = write_pairs(genes)
    copy_runs()
    verify(total)
    print(f"\npour lancer :\n  pesto --bench {TARGET.relative_to(ROOT)}/pairs.tsv")


if __name__ == "__main__":
    main()
