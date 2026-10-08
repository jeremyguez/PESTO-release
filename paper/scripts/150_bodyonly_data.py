"""Supplementary Figure 3 data: abstracts and fulltext on the body-only benchmark.

40 pairs (scripts/139_build_bench_bodyonly.py, extended by
scripts/153_extend_bench_bodyonly.py): 20 GWAS Catalog associations whose gene
is named only in the body, a table or the supplement of the GWAS article, never
in a PubMed title or abstract with the trait, and 20 negative controls. Each
pair was answered once by each arm; the result tables keep the names the arms
had when they were run (current is abstracts, currentv2-cheap is fulltext).

The verdict is the rounded mean of the distribution, as everywhere in the
paper. Writes:

  results/bodyonly_verdicts.tsv   verdict counts per arm and group
  results/bodyonly_summary.tsv    right verdicts, key PMIDs cited, cost per pair

  python3 scripts/150_bodyonly_data.py
"""
from __future__ import annotations

import csv
import statistics
import sys
from collections import Counter, defaultdict
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
score = import_module("138_score_bodyonly")

LABELS = RESULTS / "bench_bodyonly_labels.tsv"
ARMS = [("abstracts", "current"), ("fulltext", "currentv2-cheap")]
GROUP = {"gwas_body_only": "GWAS, gene in body only",
         "negative_control": "negative control"}
VERDICTS = ("Novel", "Hypothesized", "Existing", "Established")


def literature_cost(path):
    """USD per pair for the literature branch; per-pair `total` rows skipped."""
    by_pair = defaultdict(float)
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["step"] != "total" and not r["step"].startswith("ot_"):
                by_pair[(r["gene"], r["phenotype"])] += float(r["usd"])
    return list(by_pair.values())


def main():
    lab = score.labels(LABELS)
    status = {}
    with open(LABELS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            status[(r["gene"], r["phenotype"])] = GROUP[r["status"]]

    counts, summary = [], []
    for shown, stem in ARMS:
        answers = score.read_answers(RESULTS / f"bench_bodyonly_{stem}.tsv")
        rows, cited, total = score.score(lab, answers)
        c = Counter((status[(r["gene"], r["phenotype"])], r["call"]) for r in rows)
        for group in GROUP.values():
            for v in VERDICTS:
                counts.append({"arm": shown, "group": group, "verdict": v,
                               "n": c.get((group, v), 0)})
        usd = literature_cost(RESULTS / f"bench_bodyonly_{stem}.cost.tsv")
        right = Counter(status[(r["gene"], r["phenotype"])] for r in rows if r["ok"])
        size = Counter(status[(r["gene"], r["phenotype"])] for r in rows)
        summary.append({
            "arm": shown, "pairs": len(rows), "right": sum(r["ok"] for r in rows),
            "right_positive": right[GROUP["gwas_body_only"]],
            "positives": size[GROUP["gwas_body_only"]],
            "right_negative": right[GROUP["negative_control"]],
            "negatives": size[GROUP["negative_control"]],
            "pmids_cited": cited, "pmids_total": total,
            "usd_mean": f"{statistics.mean(usd):.4f}",
            "usd_sd": f"{statistics.stdev(usd):.4f}",
        })

    for path, rows in ((RESULTS / "bodyonly_verdicts.tsv", counts),
                       (RESULTS / "bodyonly_summary.tsv", summary)):
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    for s in summary:
        print(f"{s['arm']:10s} right {s['right']}/{s['pairs']} "
              f"(positives {s['right_positive']}/{s['positives']}, negatives "
              f"{s['right_negative']}/{s['negatives']}), key PMIDs cited "
              f"{s['pmids_cited']}/{s['pmids_total']}, literature "
              f"${s['usd_mean']} ± {s['usd_sd']} a pair")
    print("wrote results/bodyonly_verdicts.tsv and results/bodyonly_summary.tsv")


if __name__ == "__main__":
    main()
