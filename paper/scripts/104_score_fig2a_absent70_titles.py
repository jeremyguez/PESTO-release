"""Titles arm on the 70 extra Figure 2a-style random Absents.

Does not rerun prior or current. Writes results/fig2a_absent_extra70/titles.tsv.

Usage:
  python3 scripts/104_score_fig2a_absent70_titles.py
"""
from __future__ import annotations

import csv
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.cost import from_trace, load_table, with_total, write_table
from pesto.flow import arms, run as flow_run
from pesto.flow.types import Query

SRCS = [
    ROOT / "results" / "fig2a_absent_extra20.tsv",
    ROOT / "results" / "fig2a_absent_extra20b.tsv",
    ROOT / "results" / "fig2a_absent_extra30.tsv",
]
OUT = ROOT / "results" / "fig2a_absent_extra70"
WORKERS = 15
TITLES_FP = "1ec6267b715a"
LOCK = threading.Lock()
FIELDS = [
    "gene", "phenotype", "verdict", "lit_score",
    "found", "read", "cached", "error", "usd", "bands_usd",
]


def lit_score(dist):
    if not dist:
        return None
    return (4 * dist.get("Established", 0) + 3 * dist.get("Existing", 0)
            + 2 * dist.get("Hypothesized", 0) + 1 * dist.get("Novel", 0)) / 100


def load_tsv(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def step_usd(spent, name):
    return sum(float(r["usd"]) for r in (spent or []) if r.get("step") == name)


def pairs():
    out, seen = [], set()
    for src in SRCS:
        for r in load_tsv(src):
            k = (r["gene"], r["phenotype"])
            if k in seen:
                continue
            seen.add(k)
            out.append({"gene": r["gene"], "phenotype": r["phenotype"]})
    return out


def run_titles(todo_pairs, arm):
    path = OUT / "titles.tsv"
    cost_path = str(OUT / "titles.cost.tsv")
    runs_dir = str(OUT / "runs" / "titles")
    rows = load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in rows if not r.get("error")}
    todo = [p for p in todo_pairs if (p["gene"], p["phenotype"]) not in have]
    print(f"titles  {len(todo_pairs) - len(todo)} done, {len(todo)} to run  "
          f"[{arm.fingerprint()}]", flush=True)
    cost = load_table(cost_path)

    def work(p):
        query = Query(p["gene"], p["phenotype"])
        try:
            res = flow_run.run(arm, query, cache=True, runs_dir=runs_dir)
            dist = (res.verdict.distribution.as_dict()
                    if res.verdict.distribution else {})
            spent = with_total(from_trace(p["gene"], p["phenotype"],
                                          res.trace, res.cached),
                               p["gene"], p["phenotype"], res.cached)
            return {
                "gene": p["gene"], "phenotype": p["phenotype"],
                "verdict": res.verdict.call, "lit_score": lit_score(dist),
                "found": res.counts.get("found", ""),
                "read": res.counts.get("read", ""),
                "cached": res.cached, "error": "",
                "usd": f"{step_usd(spent, 'total'):.6f}",
                "bands_usd": f"{step_usd(spent, 'Bands'):.6f}",
            }, spent
        except Exception as exc:
            return {
                "gene": p["gene"], "phenotype": p["phenotype"],
                "verdict": "", "lit_score": "",
                "found": "", "read": "", "cached": False,
                "error": f"{type(exc).__name__}: {exc}",
                "usd": "", "bands_usd": "",
            }, None

    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, p) for p in todo]
        for fut in as_completed(futs):
            rec, spent = fut.result()
            with LOCK:
                rows.append(rec)
                if spent is not None:
                    cost[(rec["gene"], rec["phenotype"])] = spent
                save_tsv(path, rows, FIELDS)
                write_table(cost_path, cost)
                n_new += 1
                print(f"  titles {n_new}/{len(todo)}  {rec['gene']}  "
                      f"{rec.get('verdict') or rec.get('error', '')[:40]}  "
                      f"${rec.get('usd') or 0}", flush=True)
    return rows


def main():
    arm = arms.get("titles")
    if arm.fingerprint() != TITLES_FP:
        print(f"ERROR: titles fingerprint {arm.fingerprint()}", file=sys.stderr)
        return 1
    todo = pairs()
    print(f"n={len(todo)}  {WORKERS} workers  → {OUT}", flush=True)
    if len(todo) != 70:
        print(f"ERROR: expected 70 pairs, got {len(todo)}", file=sys.stderr)
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    rows = run_titles(todo, arm)
    n_err = sum(1 for r in rows if r.get("error"))
    n_ok = sum(1 for r in rows if r.get("verdict") in
               ("Novel", "Hypothesized", "Existing", "Established"))
    print(f"titles ok {n_ok}/{len(todo)}  errors {n_err}/{len(todo)}",
          flush=True)
    return 1 if n_err or n_ok < len(todo) else 0


if __name__ == "__main__":
    sys.exit(main())
