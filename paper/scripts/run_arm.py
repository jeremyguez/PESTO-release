"""Run one declared arm over a table of pairs.

What comes out carries the arm's fingerprint in every row, so two result tables
can be compared or refused on evidence rather than on memory. `--arm auto`
picks abstracts or abstracts-trait for each pair, as Figure 1 did.

  python3 scripts/run_arm.py --arm abstracts --bench ... --workers 20 --no-cache
  python3 scripts/run_arm.py --arm auto --bench results/bench_aou_brava_auto.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)

from pesto.flow import arms, run as flow_run  # noqa: E402
from pesto.flow.types import Query  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--arm", required=True,
                choices=[arms.AUTO, *sorted(arms.ARMS)])
ap.add_argument("--bench", required=True, help="TSV with gene and phenotype")
ap.add_argument("--out", help="default: results/arm_<arm>_<bench name>.json")
ap.add_argument("--workers", type=int, default=6)
ap.add_argument("--genes", nargs="+", help="only these genes")
ap.add_argument("--limit", type=int, help="only the first N pairs")
ap.add_argument("--no-cache", action="store_true",
                help="rerun even where this exact pipeline has answered the pair")
ap.add_argument("--runs-dir",
                help="where to read and write the per-pair flow cache; "
                     "an empty directory is a clean slate")
args = ap.parse_args()

auto = args.arm == arms.AUTO
arm = None if auto else arms.get(args.arm)
name = os.path.basename(args.bench).replace(".tsv", "")
out_json = args.out or os.path.join(ROOT, "results", f"arm_{args.arm}_{name}.json")
out_tsv = out_json.replace(".json", ".tsv")

with open(args.bench, encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if args.genes:
    want = {g.upper() for g in args.genes}
    rows = [r for r in rows if r["gene"].upper() in want]
if args.limit:
    rows = rows[:args.limit]


def truth_of(r):
    for key in ("gencc_class_promoted", "gencc_class", "gencc", "gencc_promoted"):
        if r.get(key):
            return r[key]
    return ""


def work(r):
    query = Query(r["gene"], r["phenotype"])
    chosen = (arms.resolve(arms.AUTO, query.phenotype) if auto else arm)
    # A benchmark's search must not move between runs, so a recorded expansion
    # is replayed where the table carries one and the arm asks for it.
    syn = [s.strip() for s in (r.get("synonyms") or "").split(";") if s.strip()]
    out = {"gene": query.gene, "phenotype": query.phenotype, "truth": truth_of(r),
           "arm": chosen.name, "fingerprint": chosen.fingerprint()}
    try:
        res = flow_run.run(chosen, query, synonyms=syn,
                           cache=not args.no_cache or bool(args.runs_dir),
                           runs_dir=args.runs_dir)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out
    v = res.verdict
    out.update({"verdict": v.call, "cached": res.cached,
                "found": res.counts.get("found", 0),
                "read": res.counts.get("read", 0),
                "justification": (v.justification or "").replace("\t", " "),
                "pmids": ",".join(v.pmids),
                "distribution": v.distribution.as_dict() if v.distribution else None,
                "detail": v.detail})
    return out


lock = threading.Lock()
done = []
started = time.time()
print(f"{len(rows)} pairs, {args.arm}"
      f"{'' if auto else ' [' + arm.fingerprint() + ']'}"
      f", {args.workers} workers")

COLUMNS = ["gene", "phenotype", "truth", "arm", "fingerprint", "verdict",
           "cached", "found", "read", "pmids", "justification", "error"]


def record(row):
    with lock:
        done.append(row)
        with open(out_json, "w", encoding="utf-8") as fh:
            json.dump(done, fh, indent=1)
        with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=COLUMNS,
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(done)
        n = len(done)
        print(f"  {n:3d}/{len(rows)}  {row['gene']:10s} "
              f"{row.get('verdict') or row.get('error', '')[:40]:14s} "
              f"{row.get('read', 0):3d} read   {time.time() - started:5.0f}s")


with ThreadPoolExecutor(max_workers=args.workers) as pool:
    for row in pool.map(work, rows):
        record(row)

by_class = {}
for row in done:
    by_class.setdefault(row["truth"], []).append(row.get("verdict"))
print()
for truth, calls in sorted(by_class.items()):
    counts = {v: calls.count(v) for v in sorted(set(calls), key=str)}
    print(f"  {truth or '(none)':12s} {len(calls):3d}  {counts}")
print(f"\nwrote {os.path.relpath(out_json, ROOT)} and its TSV")
