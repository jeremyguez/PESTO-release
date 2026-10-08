"""Extend the body-only benchmark from 16 + 6 to 20 positives and 20 controls.

The 22 pairs drawn by scripts/139_build_bench_bodyonly.py are kept as they are,
and their answers with them; redrawing from scratch would not return them,
since PubMed and Europe PMC have moved since. New pairs are drawn with the same
functions, the same per-trait generators and the same exclusions (the genes of
the earlier unpublished test set and every gene already in the benchmark):

  positives  one more for 4 of the 8 traits, chosen with a fixed seed: the
             first candidate, in that trait's seeded order, that passes the
             positive test and is not yet in the benchmark;
  controls   two for schizophrenia and blood pressure, which had none, one
             more for each of the other six traits, and one more again for 4
             traits chosen with a fixed seed, so every trait has two or three.

The original label table is kept as results/bench_bodyonly_labels_22.tsv; the
extended one replaces results/bench_bodyonly_labels.tsv and
results/bench_bodyonly.tsv, with the new pairs appended after the old ones.

  python3 scripts/153_extend_bench_bodyonly.py
"""
from __future__ import annotations

import csv
import random
import shutil
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from importlib import import_module
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
build = import_module("139_build_bench_bodyonly")

ROOT = build.ROOT
LABELS = build.OUT
BENCH = build.BENCH
KEPT = ROOT / "results" / "bench_bodyonly_labels_22.tsv"
NEW = ROOT / "results" / "bench_bodyonly_new18.tsv"
TARGET = 20
COLS = ["gene", "phenotype", "type", "accept", "key_pmids", "status", "note", "split"]


def main():
    with LABELS.open(encoding="utf-8") as fh:
        old = list(csv.DictReader(fh, delimiter="\t"))
    if len(old) != 22:
        raise SystemExit(f"expected the 22 published pairs, found {len(old)}")
    build.CAT.update(build.load_catalog())
    names = build.hgnc()
    used = build.already_used() | {r["gene"] for r in old}
    by_trait = {t[1]: t for t in build.TRAITS}

    # Positives: one more for four traits.
    n_pos = sum(r["status"] == "gwas_body_only" for r in old)
    pos_traits = random.Random(f"{build.SEED}:extend:positive").sample(
        build.TRAITS, TARGET - n_pos)
    build.POSITIVES_PER_TRAIT = 1
    with ThreadPoolExecutor(len(pos_traits)) as ex:
        parts = list(ex.map(lambda t: build.positives(*t, names, used), pos_traits))
    new = [r for part in parts for r in part]
    used |= {r["gene"] for r in new}

    # Controls: how many each trait still needs.
    have = Counter(r["phenotype"] for r in old if r["status"] == "negative_control")
    want = {t[1]: 2 for t in build.TRAITS}
    for t in random.Random(f"{build.SEED}:extend:negative").sample(
            sorted(want), TARGET - 2 * len(want)):
        want[t] += 1
    need = {p: want[p] - have.get(p, 0) for p in want}

    def controls(pheno):
        rows = []
        for _ in range(need[pheno]):
            got = build.negative(*by_trait[pheno], names, used | {r["gene"] for r in rows})
            rows += got
        return rows

    with ThreadPoolExecutor(len(need)) as ex:
        for part in ex.map(controls, sorted(need)):
            for r in part:
                if r["gene"] not in used:
                    used.add(r["gene"])
                    new.append(r)

    shutil.copyfile(LABELS, KEPT)
    rows = old + new
    for path, cols, data in ((LABELS, COLS, rows),
                             (NEW, ["gene", "phenotype"], new),
                             (BENCH, ["gene", "phenotype"], rows)):
        with path.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(data)
    c = Counter(r["status"] for r in rows)
    print(f"{len(new)} new pairs; benchmark now {c['gwas_body_only']} positives, "
          f"{c['negative_control']} controls -> {LABELS.relative_to(ROOT)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
