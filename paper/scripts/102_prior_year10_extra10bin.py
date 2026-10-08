"""Ten more Limited pairs per 4-year bin on top of prior_year10.

Skips genes already in results/prior_year10/sample.tsv. Seed 2030. Scores
prior then current and appends to sample.tsv, knowledge.tsv, pesto.json
and pesto.tsv only after both arms finish.

Usage:
  python3 scripts/102_prior_year10_extra10bin.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import random
import sys
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.flow import arms, run as flow_run
from pesto.flow.types import Query

spec = importlib.util.spec_from_file_location(
    "s80", ROOT / "scripts" / "80_prior_year_sample.py")
s80 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s80)

OUT = ROOT / "results" / "prior_year10"
EXTRA = OUT / "extra10bin.tsv"
N_PER_BIN = 10
SEED = 2030
WORKERS = 20
PRIOR_FP = "2b18abb989c5"
CURRENT_FP = "ff136250f372"
LOCK = threading.Lock()
BINS = (
    ("2010-2013", range(2010, 2014)),
    ("2014-2017", range(2014, 2018)),
    ("2018-2021", range(2018, 2022)),
    ("2022-2025", range(2022, 2026)),
)
SAMPLE_FIELDS = [
    "year", "gene", "phenotype", "disease_curie", "n_pmids", "oldest_pmid",
]
PESTO_TSV_FIELDS = [
    "gene", "phenotype", "truth", "arm", "fingerprint", "verdict",
    "cached", "found", "read", "pmids", "justification", "error",
]


def bin_of(year):
    y = int(year)
    if y <= 2013:
        return "2010-2013"
    if y <= 2017:
        return "2014-2017"
    if y <= 2021:
        return "2018-2021"
    return "2022-2025"


def draw(dated, spent_keys, spent_genes):
    by_bin = defaultdict(list)
    for p in dated:
        if not (s80.YEAR0 <= p["year"] <= s80.YEAR1):
            continue
        if (p["gene"], p["phenotype"]) in spent_keys or p["gene"] in spent_genes:
            continue
        by_bin[bin_of(p["year"])].append(p)
    rng = random.Random(SEED)
    sample = []
    print("pool / drawn by bin:", flush=True)
    for label, _years in BINS:
        pool = list(by_bin.get(label, []))
        rng.shuffle(pool)
        take = pool[:N_PER_BIN]
        print(f"  {label}: {len(pool)} -> {len(take)}", flush=True)
        if len(take) < N_PER_BIN:
            sys.exit(f"only {len(take)} left in {label}")
        sample.extend(take)
    return sample


def sample_row(p):
    return {k: p[k] for k in SAMPLE_FIELDS}


def run_prior(todo):
    path = OUT / "extra10bin_knowledge.tsv"
    rows = s80.load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in rows
            if r.get("verdict") in ("Novel", "Hypothesized", "Existing",
                                    "Established")}
    need = [p for p in todo if (p["gene"], p["phenotype"]) not in have]
    print(f"prior {len(todo) - len(need)} done, {len(need)} to run", flush=True)
    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(s80.knowledge_one, p["gene"], p["phenotype"]): p
                for p in need}
        for fut in as_completed(futs):
            p = futs[fut]
            try:
                one = fut.result()
            except Exception as exc:
                one = {"verdict": "", "lit_score": "", "error": str(exc)}
                print(f"  FAIL {p['gene']}: {exc}", flush=True)
            rec = {**sample_row(p), **one}
            rec.pop("pmids", None)
            with LOCK:
                rows.append(rec)
                s80.save_tsv(path, rows, s80.KNOW_FIELDS)
                n_new += 1
                usd = rec.get("usd") or 0
                print(f"  prior {n_new}/{len(need)}  {p['year']}  {p['gene']}  "
                      f"{rec.get('verdict')}  ${float(usd or 0):.3f}",
                      flush=True)
    return rows


def run_current(todo):
    path = OUT / "extra10bin_pesto.json"
    tsv_path = OUT / "extra10bin_pesto.tsv"
    if path.exists():
        rows = json.loads(path.read_text())
    else:
        rows = []
    have = {(r["gene"], r["phenotype"]) for r in rows
            if r.get("verdict") and not r.get("error")}
    need = [p for p in todo if (p["gene"], p["phenotype"]) not in have]
    arm = arms.get("current")
    runs_dir = str(OUT / "runs" / "flow")
    print(f"current {len(todo) - len(need)} done, {len(need)} to run  "
          f"[{arm.fingerprint()}]", flush=True)

    def work(p):
        query = Query(p["gene"], p["phenotype"])
        try:
            res = flow_run.run(arm, query, cache=True, runs_dir=runs_dir)
            v = res.verdict
            dist = v.distribution.as_dict() if v.distribution else None
            return {
                "gene": p["gene"], "phenotype": p["phenotype"], "truth": "",
                "arm": "current", "fingerprint": CURRENT_FP,
                "verdict": v.call, "cached": res.cached,
                "found": res.counts.get("found", 0),
                "read": res.counts.get("read", 0),
                "justification": (v.justification or "").replace("\t", " "),
                "pmids": ",".join(v.pmids),
                "distribution": dist, "detail": v.detail, "error": "",
            }
        except Exception as exc:
            return {
                "gene": p["gene"], "phenotype": p["phenotype"], "truth": "",
                "arm": "current", "fingerprint": CURRENT_FP,
                "verdict": "", "cached": False, "found": "", "read": "",
                "justification": "", "pmids": "", "distribution": None,
                "detail": None, "error": f"{type(exc).__name__}: {exc}",
            }

    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, p) for p in need]
        for fut in as_completed(futs):
            rec = fut.result()
            with LOCK:
                rows.append(rec)
                path.write_text(json.dumps(rows, indent=1), encoding="utf-8")
                with tsv_path.open("w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, PESTO_TSV_FIELDS, delimiter="\t",
                                       extrasaction="ignore")
                    w.writeheader()
                    w.writerows(rows)
                n_new += 1
                print(f"  current {n_new}/{len(need)}  {rec['gene']}  "
                      f"{rec.get('verdict') or rec.get('error', '')[:40]}",
                      flush=True)
    return rows


def append_tsv(path, rows, fields):
    existing = s80.load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in existing}
    add = [r for r in rows if (r["gene"], r["phenotype"]) not in have]
    s80.save_tsv(path, existing + add, fields)
    return len(add)


def append_json(path, rows):
    existing = json.loads(path.read_text()) if path.exists() else []
    have = {(r["gene"], r["phenotype"]) for r in existing}
    add = [r for r in rows if (r["gene"], r["phenotype"]) not in have]
    existing.extend(add)
    path.write_text(json.dumps(existing, indent=1), encoding="utf-8")
    return len(add)


def main():
    prior, current = arms.get("prior"), arms.get("current")
    if prior.fingerprint() != PRIOR_FP:
        print(f"ERROR: prior fingerprint {prior.fingerprint()}", file=sys.stderr)
        return 1
    if current.fingerprint() != CURRENT_FP:
        print(f"ERROR: current fingerprint {current.fingerprint()}",
              file=sys.stderr)
        return 1

    have_sample = s80.load_tsv(OUT / "sample.tsv")
    spent_keys = {(r["gene"], r["phenotype"]) for r in have_sample}
    spent_genes = {r["gene"] for r in have_sample}
    pairs = s80.limited_pairs()
    dated = s80.date_pairs(pairs, s80.load_pmid_years())
    extra = draw(dated, spent_keys, spent_genes)
    s80.save_tsv(EXTRA, [sample_row(p) for p in extra], SAMPLE_FIELDS)
    print(f"wrote {EXTRA} n={len(extra)}", flush=True)

    know = run_prior(extra)
    pesto = run_current(extra)
    n_err = sum(1 for r in pesto if r.get("error"))
    n_k = sum(1 for r in know if r.get("verdict") in
              ("Novel", "Hypothesized", "Existing", "Established"))
    if n_err or n_k < len(extra):
        print(f"prior ok {n_k}/{len(extra)}  pesto errors {n_err}/{len(extra)} "
              f"— not merging", flush=True)
        return 1

    n_s = append_tsv(OUT / "sample.tsv", [sample_row(p) for p in extra],
                     SAMPLE_FIELDS)
    n_k2 = append_tsv(OUT / "knowledge.tsv", know, s80.KNOW_FIELDS)
    n_pj = append_json(OUT / "pesto.json", pesto)
    n_pt = append_tsv(OUT / "pesto.tsv", pesto, PESTO_TSV_FIELDS)
    print(f"merged +{n_s} sample  +{n_k2} knowledge  "
          f"+{n_pj} pesto.json  +{n_pt} pesto.tsv", flush=True)
    for label, _years in BINS:
        xs = [p for p in extra if bin_of(p["year"]) == label]
        print(f"  {label}: " + ", ".join(p["gene"] for p in xs), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
