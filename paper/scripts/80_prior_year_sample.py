"""Closed-book prior on 10 Limited pairs a year, 2010–2025.

Oldest cited PMID year, reappraisal PMIDs ignored, the nine HCM G2P
reappraisals dropped. Writes results/prior_year10/. Prior only, 20 workers.
2024 takes 17 (ten plus seven more) so the 2022–2025 bin matches n=40.

Usage:
  python3 scripts/80_prior_year_sample.py
"""
from __future__ import annotations

import csv
import os
import random
import re
import sys
import threading
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.cost import capturing
from pesto.flow import arms
from pesto.flow.steps import Argmax
from pesto.flow.types import Query
from pesto.harness import MODELS
from pesto.services.llm_service import call_llm_with_usage
from pesto.services.pubmed_service import (
    BASE_URL_NCBI, NCBI_API_KEY, make_api_request_with_retry,
)
from pesto.utils.helpers import load_prompt

OUT = ROOT / "results" / "prior_year10"
PMID_YEARS = ROOT / "results" / "pmid_pubyears.tsv"
BENCH = OUT / "sample.tsv"
KNOW = OUT / "knowledge.tsv"
WORKERS = 20
SEED = 2026
PER_YEAR = 10
EXTRA_2024 = 7
YEAR0, YEAR1 = 2010, 2025
PROMPT = load_prompt("knowledge_probabilities")
LOCK = threading.Lock()
PMID_RE = re.compile(r"(\d{5,9})")
REAPP = {"39971408", "33831308", "31983240", "29959160"}
HCM = {
    "KLF10", "NEXN", "OBSCN", "PDLIM3", "RBM20",
    "RPS6KB1", "RYR2", "TMPO", "TTN",
}
PRIOR_FP = "2b18abb989c5"
KNOW_FIELDS = [
    "year", "gene", "phenotype", "disease_curie", "n_pmids", "oldest_pmid",
    "verdict", "lit_score",
    "p_established", "p_existing", "p_hypothesized", "p_novel",
    "justification", "input_tokens", "output_tokens", "usd", "calls",
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


def load_pmid_years():
    years = {}
    for r in load_tsv(PMID_YEARS):
        p, y = r.get("pmid"), r.get("year")
        if p and (y or "").isdigit():
            years[p] = int(y)
    return years


def save_pmid_years(years):
    rows = [{"pmid": p, "year": y} for p, y in sorted(years.items())]
    save_tsv(PMID_YEARS, rows, ["pmid", "year"])


def fetch_years(pmids):
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    resp = make_api_request_with_retry(
        f"{BASE_URL_NCBI}esummary.fcgi", params, timeout=60)
    if resp is None:
        return {}
    out = {}
    root = ET.fromstring(resp.content)
    for doc in root.findall(".//DocSum"):
        pid_el = doc.find("Id")
        pid = pid_el.text.strip() if pid_el is not None and pid_el.text else None
        pubdate = ""
        for item in doc.findall("Item"):
            if item.get("Name") == "PubDate":
                pubdate = (item.text or "").strip()
        if pid and pubdate[:4].isdigit():
            out[pid] = int(pubdate[:4])
    return out


def limited_pairs():
    subs = defaultdict(list)
    with (ROOT / "data/raw/gencc-submissions.tsv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = (r.get("gene_symbol") or "").strip()
            curie = (r.get("disease_curie") or "").strip()
            if gene and curie:
                subs[(gene, curie)].append(r)
    rank = {"Definitive": 0, "Strong": 1, "Moderate": 2, "Limited": 3}
    out = []
    for (gene, curie), lst in subs.items():
        best = None
        for r in lst:
            c = r.get("classification_title")
            if c in rank and (best is None or rank[c] < rank[best]):
                best = c
        if best != "Limited":
            continue
        titles = Counter()
        pmids = []
        for r in lst:
            if r.get("classification_title") != "Limited":
                continue
            t = (r.get("disease_title") or "").strip()
            if t:
                titles[t] += 1
            for p in PMID_RE.findall(r.get("submitted_as_pmids") or ""):
                if p not in pmids:
                    pmids.append(p)
        phenotype = titles.most_common(1)[0][0] if titles else ""
        if gene in HCM and "hypertrophic cardiomyopathy" in phenotype.lower():
            continue
        usable = [p for p in pmids if p not in REAPP]
        if not usable:
            continue
        out.append({
            "gene": gene, "phenotype": phenotype, "disease_curie": curie,
            "pmids": usable,
        })
    return out


def date_pairs(pairs, years):
    need = []
    seen = set()
    for p in pairs:
        for pmid in p["pmids"]:
            if pmid not in years and pmid not in seen:
                seen.add(pmid)
                need.append(pmid)
    print(f"PMID years cached {len(years)}, to fetch {len(need)}", flush=True)
    for i in range(0, len(need), 200):
        chunk = need[i:i + 200]
        got = fetch_years(chunk)
        years.update(got)
        print(f"  batch {i // 200 + 1} +{len(got)}", flush=True)
        save_pmid_years(years)
    dated = []
    for p in pairs:
        ys = [(pmid, years[pmid]) for pmid in p["pmids"] if pmid in years]
        if not ys:
            continue
        pmid, y = min(ys, key=lambda t: t[1])
        dated.append({**p, "year": y, "oldest_pmid": pmid, "n_pmids": len(p["pmids"])})
    return dated


def draw(dated):
    by_year = defaultdict(list)
    for p in dated:
        if YEAR0 <= p["year"] <= YEAR1:
            by_year[p["year"]].append(p)
    rng = random.Random(SEED)
    sample = []
    print("pool / drawn:", flush=True)
    for y in range(YEAR0, YEAR1 + 1):
        pool = by_year.get(y, [])
        rng.shuffle(pool)
        n = PER_YEAR + (EXTRA_2024 if y == 2024 else 0)
        take = pool[:n]
        print(f"  {y}: {len(pool)} -> {len(take)}", flush=True)
        sample.extend(take)
    return sample


def knowledge_one(gene, phenotype):
    text = PROMPT.format(gene_name=gene, phenotype=phenotype)
    with capturing() as cap:
        resp = call_llm_with_usage(
            text, MODELS["opus5"], 0.0, agent_name="knowledge_prior") or {}
    raw = resp.get("text", "") or ""
    verdict = Argmax().run(Query(gene, phenotype), raw)
    spent = cap.summary(MODELS["opus5"])
    dist = verdict.distribution.as_dict() if verdict.distribution else {}
    return {
        "verdict": verdict.call, "lit_score": lit_score(dist),
        "p_established": dist.get("Established", ""),
        "p_existing": dist.get("Existing", ""),
        "p_hypothesized": dist.get("Hypothesized", ""),
        "p_novel": dist.get("Novel", ""),
        "justification": (verdict.justification or "").replace("\t", " "),
        **spent,
    }


def run_prior(sample):
    have = {(r["gene"], r["phenotype"]): r for r in load_tsv(KNOW)
            if r.get("verdict") in ("Novel", "Hypothesized", "Existing", "Established")}
    rows = list(have.values())
    todo = [p for p in sample if (p["gene"], p["phenotype"]) not in have]
    print(f"prior {len(sample) - len(todo)} cached, {len(todo)} to run",
          flush=True)
    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(knowledge_one, p["gene"], p["phenotype"]): p
                for p in todo}
        for fut in as_completed(futs):
            p = futs[fut]
            try:
                one = fut.result()
            except Exception as exc:
                one = {"verdict": "", "lit_score": "", "error": str(exc)}
                print(f"  FAIL {p['gene']}: {exc}", flush=True)
            rec = {**p, **one}
            rec.pop("pmids", None)
            with LOCK:
                rows.append(rec)
                save_tsv(KNOW, rows, KNOW_FIELDS)
                n_new += 1
                usd = rec.get("usd") or 0
                print(f"  prior {n_new}/{len(todo)}  {p['year']}  {p['gene']}  "
                      f"{rec.get('verdict')}  ${float(usd or 0):.3f}",
                      flush=True)
    return rows


def report(rows):
    ee = lambda v: v in ("Existing", "Established")
    print("\n=== prior by oldest-cited year ===", flush=True)
    print("year\tn\tEE\tpct\tmean", flush=True)
    by = defaultdict(list)
    for r in rows:
        try:
            y = int(r["year"])
        except (TypeError, ValueError, KeyError):
            continue
        if r.get("verdict") in ("Novel", "Hypothesized", "Existing", "Established"):
            by[y].append(r)
    for y in range(YEAR0, YEAR1 + 1):
        xs = by.get(y, [])
        if not xs:
            print(f"{y}\t0", flush=True)
            continue
        n_ee = sum(ee(r["verdict"]) for r in xs)
        scores = [float(r["lit_score"]) for r in xs if r.get("lit_score") not in ("", None)]
        mean = sum(scores) / len(scores) if scores else float("nan")
        print(f"{y}\t{len(xs)}\t{n_ee}/{len(xs)}\t{100 * n_ee / len(xs):.0f}%\t{mean:.2f}",
              flush=True)
    old = [r for y, xs in by.items() if y < 2023 for r in xs]
    new = by.get(2023, []) + by.get(2024, []) + by.get(2025, [])
    for label, xs in (("<2023", old), ("2023-2025", new)):
        if not xs:
            continue
        n_ee = sum(ee(r["verdict"]) for r in xs)
        scores = [float(r["lit_score"]) for r in xs if r.get("lit_score") not in ("", None)]
        print(f"{label}\t{len(xs)}\t{n_ee}/{len(xs)}\t{100 * n_ee / len(xs):.0f}%\t"
              f"{sum(scores) / len(scores):.2f}", flush=True)


def main():
    prior = arms.get("prior")
    if prior.fingerprint() != PRIOR_FP:
        print(f"ERROR: prior fingerprint {prior.fingerprint()} != {PRIOR_FP}",
              file=sys.stderr)
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = limited_pairs()
    print(f"Limited pairs after HCM/reappraisal drop: {len(pairs)}", flush=True)
    dated = date_pairs(pairs, load_pmid_years())
    sample = draw(dated)
    save_tsv(BENCH, [{k: p[k] for k in
                      ("year", "gene", "phenotype", "disease_curie",
                       "n_pmids", "oldest_pmid")} for p in sample],
             ["year", "gene", "phenotype", "disease_curie", "n_pmids", "oldest_pmid"])
    print(f"wrote {BENCH} n={len(sample)}", flush=True)
    rows = run_prior(sample)
    report(rows)
    usd = sum(float(r["usd"]) for r in rows if r.get("usd") not in ("", None))
    print(f"total ${usd:.2f}  wrote {KNOW}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
