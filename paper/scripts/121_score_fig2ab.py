"""Novelty-score the Figure 2a and 2b pairs, reusing `current` articles.

Seeds scores already computed on the 2023–2025 AP Absent set.

  PESTO_PROJECT_ROOT=$PWD python3 scripts/121_score_fig2ab.py
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.flow import arms
from pesto.flow.steps import NoveltyScoreRead
from pesto.flow.store import slug
from pesto.flow.types import Corpus, Query, FIELDS as ARTICLE_FIELDS

CURRENT_FP = "ff136250f372"
WORKERS = 20
BENCH = ROOT / "results" / "bench_fig2ab_score.tsv"
OUT_JSON = ROOT / "results" / "arm_current-score_fig2ab.json"
OUT_TSV = ROOT / "results" / "arm_current-score_fig2ab.tsv"
SEED = ROOT / "results" / "arm_current-score_fig2d_ap_2023_2025.json"
LOCK = threading.Lock()
COLUMNS = ["gene", "phenotype", "panel", "group", "novelty", "evidence",
           "read", "error"]
OK = {"Novel", "Hypothesized", "Existing", "Established"}
ASD_GROUP = {
    "autism known": "ASD known",
    "ID known": "ASD candidate",
    "autism candidate": "ASD candidate",
    "novel candidate": "novel candidate",
}
NOVELTY = re.compile(
    r"<novelty>\s*([0-9]+(?:\.[0-9]+)?)\s*</novelty>", re.I)
FLOW_INDEX = None


def load_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_novelty(text):
    m = NOVELTY.search(text or "")
    if not m:
        return None
    value = float(m.group(1))
    if value < 0 or value > 100:
        return None
    return value


def flow_index():
    global FLOW_INDEX
    if FLOW_INDEX is None:
        FLOW_INDEX = {}
        for folder in ROOT.glob("**/runs/flow"):
            for path in folder.glob(f"*__{CURRENT_FP}.json"):
                FLOW_INDEX[path.name] = path
    return FLOW_INDEX


def current_run(query):
    name = f"{slug(query)}__{CURRENT_FP}.json"
    path = flow_index().get(name)
    if path is None:
        raise FileNotFoundError(name)
    return path


def pairs():
    rows, seen = [], set()

    def add(panel, group, gene, phenotype):
        key = (gene.upper(), (phenotype or "").strip().lower())
        if key in seen:
            return
        seen.add(key)
        rows.append({"gene": gene, "phenotype": phenotype,
                     "panel": panel, "group": group})

    for r in load_tsv(ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv"):
        if r.get("gencc_class") in {
                "Definitive", "Strong", "Moderate", "Limited", "Absent"}:
            add("2a", r["gencc_class"], r["gene"], r["phenotype"])
    for folder in ("fig2a_absent_extra20", "fig2a_absent_extra20b",
                   "fig2a_absent_extra30"):
        for r in load_tsv(ROOT / "results" / folder / "pesto.tsv"):
            if not r.get("error"):
                add("2a", "Absent", r["gene"], r["phenotype"])

    asd = load_tsv(ROOT / "results" / "bench_asd_current.tsv")
    have = set()
    for r in asd:
        group = ASD_GROUP.get(r.get("classification"))
        if not group:
            continue
        add("2b", group, r["gene"], r["phenotype"])
        have.add(r["gene"].upper())
    extra101 = [r for r in load_tsv(ROOT / "benchmark" / "asd_extra101" / "pesto.tsv")
                if r.get("set") == "random" and r.get("lit_mean") in OK]
    extra101.sort(key=lambda r: r["gene"].upper())
    for r in extra101[:50]:
        if r["gene"].upper() in have:
            continue
        add("2b", "not in ASC", r["gene"], r["phenotype"])
        have.add(r["gene"].upper())
    for r in load_tsv(ROOT / "benchmark" / "extra50_random" / "pesto.tsv"):
        if r.get("phenotype") != "Autism Spectrum Disorder" or r.get("lit_mean") not in OK:
            continue
        if r["gene"].upper() in have:
            continue
        add("2b", "not in ASC", r["gene"], r["phenotype"])
        have.add(r["gene"].upper())
    return rows


def save_progress(done):
    with OUT_JSON.open("w", encoding="utf-8") as fh:
        json.dump(done, fh, indent=1)
    with OUT_TSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(done)


def work(row, reader):
    query = Query(row["gene"], row["phenotype"])
    out = {"gene": query.gene, "phenotype": query.phenotype,
           "panel": row["panel"], "group": row["group"]}
    try:
        raw = json.loads(current_run(query).read_text())
        corpus = Corpus.of(raw.get("articles") or [], ARTICLE_FIELDS)
        text, _ = reader.run(query, corpus)
        novelty = parse_novelty(text)
        if novelty is None:
            text, _ = reader.run(query, corpus)
            novelty = parse_novelty(text)
        if novelty is None:
            out["error"] = "PARSE_FAIL"
            out["raw"] = text
            out["read"] = len(corpus)
            return out
        out.update({
            "novelty": novelty,
            "evidence": 100.0 - novelty,
            "read": len(corpus),
            "raw": text,
        })
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def main():
    current = arms.get("current")
    if current.fingerprint() != CURRENT_FP:
        print(f"ERROR: current fingerprint moved to {current.fingerprint()}",
              file=sys.stderr)
        return 1
    rows = pairs()
    n2a = sum(1 for r in rows if r["panel"] == "2a")
    n2b = sum(1 for r in rows if r["panel"] == "2b")
    print("groups", dict(Counter((r["panel"], r["group"]) for r in rows)))
    if n2a != 220 or n2b != 353:
        print(f"ERROR: expected 220+353, got {n2a}+{n2b}", file=sys.stderr)
        return 1
    BENCH.parent.mkdir(parents=True, exist_ok=True)
    with BENCH.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "panel", "group"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    flow_index()
    for r in rows:
        current_run(Query(r["gene"], r["phenotype"]))

    seed = {}
    if SEED.exists():
        for r in json.loads(SEED.read_text()):
            if r.get("novelty") is None or r.get("error"):
                continue
            seed[(r["gene"].upper(), r["phenotype"].strip().lower())] = r

    have = []
    if OUT_JSON.exists():
        have = json.loads(OUT_JSON.read_text())
    done_keys = {(r["gene"].upper(), r["phenotype"].strip().lower())
                 for r in have if r.get("novelty") is not None
                 and not r.get("error")}
    seeded = []
    todo = []
    for r in rows:
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        if key in done_keys:
            continue
        if key in seed:
            s = seed[key]
            seeded.append({
                "gene": r["gene"], "phenotype": r["phenotype"],
                "panel": r["panel"], "group": r["group"],
                "novelty": s["novelty"], "evidence": s.get("evidence",
                    100.0 - float(s["novelty"])),
                "read": s.get("read", ""), "cached": True,
            })
            continue
        todo.append(r)

    done = list(have) + seeded
    if seeded:
        save_progress(done)
    print(f"{len(rows)} pairs, {len(done_keys)} cached, {len(seeded)} seeded, "
          f"{len(todo)} to run, {WORKERS} workers", flush=True)
    if not todo:
        save_progress(done)
        print(f"wrote {OUT_TSV.relative_to(ROOT)}")
        return 0

    reader = NoveltyScoreRead("opus5")
    started = time.time()
    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, r, reader) for r in todo]
        for fut in as_completed(futs):
            rec = fut.result()
            with LOCK:
                done.append(rec)
                save_progress(done)
                n_new += 1
                shown = rec.get("error") or f"n={rec.get('novelty')}"
                print(f"  {n_new:3d}/{len(todo)}  {rec['gene']:10s} "
                      f"{shown:16s}  {time.time() - started:5.0f}s",
                      flush=True)
    n_ok = sum(1 for r in done if r.get("novelty") is not None)
    n_bad = sum(1 for r in done if r.get("error"))
    print(f"\nparsed {n_ok}/{len(done)}  errors {n_bad}")
    print(f"wrote {OUT_JSON.relative_to(ROOT)} and its TSV")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
