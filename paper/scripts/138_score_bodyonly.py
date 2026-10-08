"""Score an arm against the body-only benchmark.

Reads results/bench_bodyonly_labels.tsv for the labels and a pesto --bench
output table, and reports three things, because a verdict alone hides the
failure this benchmark was built for:

  verdict     the call is in the accepted set for the pair
  evidence    the PMIDs a pair's label rests on are cited in the answer
  adjacency   how far a wrong call is, on the Novel(1)..Established(4) ladder

The verdict is the rounded mean of the 100-point distribution (`lit_mean`), the
rule used throughout the paper. Evidence recall is the number that moves when
retrieval improves: an arm can land on a plausible verdict without having seen
the article the label rests on. Pairs whose label carries no PMIDs (the
negative controls) are skipped in the evidence column rather than counted as
misses.

  python3 scripts/138_score_bodyonly.py results/bench_bodyonly_current.tsv
  python3 scripts/138_score_bodyonly.py a.tsv b.tsv --names abstracts fulltext
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "results" / "bench_bodyonly_labels.tsv"
LADDER = {"Novel": 1, "Hypothesized": 2, "Existing": 3, "Established": 4}


def labels(path=SPLIT):
    out = {}
    for r in csv.DictReader(Path(path).open(encoding="utf-8"), delimiter="\t"):
        out[(r["gene"], r["phenotype"])] = {
            "accept": set(r["accept"].split("|")),
            "pmids": {p for p in r["key_pmids"].split(";") if p.strip()},
            "type": r["type"], "split": r["split"], "note": r["note"],
        }
    return out


def read_answers(path):
    out = {}
    for r in csv.DictReader(Path(path).open(encoding="utf-8"), delimiter="\t"):
        cited = set(re.findall(r"\b\d{7,8}\b", r.get("justification", "")))
        out[(r["gene"], r["phenotype"])] = {
            "call": r.get("lit_mean", ""),
            "cited": cited,
            "found": int(r.get("found") or 0),
            "read": int(r.get("read") or 0),
            "arm": r.get("arm", ""),
            "justification": r.get("justification", ""),
        }
    return out


def score(lab, ans):
    rows, ev_hit, ev_tot = [], 0, 0
    for pair, meta in lab.items():
        a = ans.get(pair)
        if a is None:
            continue
        ok = a["call"] in meta["accept"]
        dist = 0 if ok else min(abs(LADDER[a["call"]] - LADDER[t])
                                for t in meta["accept"]) if a["call"] in LADDER else 9
        missing = meta["pmids"] - a["cited"]
        if meta["pmids"]:
            ev_tot += len(meta["pmids"])
            ev_hit += len(meta["pmids"] & a["cited"])
        rows.append({"gene": pair[0], "phenotype": pair[1], "type": meta["type"],
                     "accept": "|".join(sorted(meta["accept"])), "call": a["call"],
                     "ok": ok, "dist": dist, "missing": ";".join(sorted(missing)),
                     "found": a["found"], "read": a["read"], "arm": a["arm"]})
    return rows, ev_hit, ev_tot


def report(name, rows, ev_hit, ev_tot):
    n = len(rows)
    ok = sum(r["ok"] for r in rows)
    print(f"\n=== {name}: {ok}/{n} verdicts in the accepted set "
          f"({100 * ok / n:.0f}%), evidence {ev_hit}/{ev_tot} cited")
    for t in ("disease", "trait", "selection"):
        sub = [r for r in rows if r["type"] == t]
        if sub:
            print(f"    {t:10s} {sum(r['ok'] for r in sub)}/{len(sub)}")
    bad = [r for r in rows if not r["ok"]]
    if bad:
        print("  misses:")
        for r in sorted(bad, key=lambda r: -r["dist"]):
            print(f"    {r['gene']:9s} {r['phenotype'][:32]:32s} want {r['accept']:22s} "
                  f"got {r['call']:13s} d={r['dist']}")
    gaps = [r for r in rows if r["missing"]]
    if gaps:
        print("  evidence not cited:")
        for r in gaps:
            print(f"    {r['gene']:9s} {r['phenotype'][:32]:32s} {r['missing']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tables", nargs="+")
    ap.add_argument("--names", nargs="*", default=None)
    ap.add_argument("--labels", default=SPLIT,
                    help="label table (default: results/bench_bodyonly_labels.tsv)")
    args = ap.parse_args()
    lab = labels(args.labels)
    names = args.names or [Path(t).stem for t in args.tables]
    for name, table in zip(names, args.tables):
        rows, hit, tot = score(lab, read_answers(table))
        report(name, rows, hit, tot)


if __name__ == "__main__":
    main()
