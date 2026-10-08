"""Fill each 3-year bin (2011-2013 … 2023-2025) to 40 Limited pairs.

Leaves 2010 aside. Skips genes already in results/prior_year10/sample.tsv.
Seed 2031. Scores prior then current and appends only after both finish.

Usage:
  python3 scripts/103_prior_year10_fill3year.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import random
import sys
import threading
from collections import Counter, defaultdict
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
EXTRA = OUT / "extra3year.tsv"
TARGET = 40
SEED = 2031
WORKERS = 20
PRIOR_FP = "2b18abb989c5"
CURRENT_FP = "ff136250f372"
LOCK = threading.Lock()
BINS = (
    ("2011-2013", range(2011, 2014)),
    ("2014-2016", range(2014, 2017)),
    ("2017-2019", range(2017, 2020)),
    ("2020-2022", range(2020, 2023)),
    ("2023-2025", range(2023, 2026)),
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
    if y <= 2010:
        return None
    if y <= 2013:
        return "2011-2013"
    if y <= 2016:
        return "2014-2016"
    if y <= 2019:
        return "2017-2019"
    if y <= 2022:
        return "2020-2022"
    return "2023-2025"


def draw(dated, have_sample):
    spent_keys = {(r["gene"], r["phenotype"]) for r in have_sample}
    spent_genes = {r["gene"] for r in have_sample}
    have_bin = Counter(bin_of(r["year"]) for r in have_sample if bin_of(r["year"]))
    by_bin = defaultdict(list)
    for p in dated:
        b = bin_of(p["year"])
        if not b:
            continue
        if (p["gene"], p["phenotype"]) in spent_keys or p["gene"] in spent_genes:
            continue
        by_bin[b].append(p)
    rng = random.Random(SEED)
    sample, taken_genes = [], set()
    print("have / need / pool / drawn by bin:", flush=True)
    for label, _years in BINS:
        need = TARGET - have_bin.get(label, 0)
        pool = [p for p in by_bin.get(label, []) if p["gene"] not in taken_genes]
        rng.shuffle(pool)
        take = pool[:max(need, 0)]
        print(f"  {label}: have {have_bin.get(label, 0)} need {need} "
              f"pool {len(pool)} -> {len(take)}", flush=True)
        if need > 0 and len(take) < need:
            sys.exit(f"only {len(take)} left in {label}")
        for p in take:
            taken_genes.add(p["gene"])
        sample.extend(take)
    return sample


def sample_row(p):
    return {k: p[k] for k in SAMPLE_FIELDS}


def run_prior(todo):
    path = OUT / "extra3year_knowledge.tsv"
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
    path = OUT / "extra3year_pesto.json"
    tsv_path = OUT / "extra3year_pesto.tsv"
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
    pairs = s80.limited_pairs()
    dated = s80.date_pairs(pairs, s80.load_pmid_years())
    extra = draw(dated, have_sample)
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
