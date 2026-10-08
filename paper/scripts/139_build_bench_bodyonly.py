"""Build a held-out benchmark of pairs whose evidence sits only in article bodies.

A test set of well-known pairs cannot separate fulltext from abstracts: none of
them depends on a gene named only in a table or a results section, which is the
one thing fulltext adds. This builds pairs that do, without running any pesto
arm, so the set is not shaped by the pipeline it judges.

Positives. A gene mapped to a genome-wide significant association (p <= 5e-8)
in the GWAS Catalog, with the gene alone at that locus, such that
  - PubMed finds no title or abstract naming the gene (or an alias) together
    with the trait, so the abstract route cannot reach it, and
  - Europe PMC indexes the gene in the body, a table or the supplement of the
    GWAS article itself, so the evidence exists in text a reader could find.
Accepted verdicts Existing|Established: one credible human association study
at least. The key PMIDs are the catalog articles that name the gene.

Negatives. A protein-coding gene with no catalog association to the trait, no
PubMed title or abstract pairing, and no Europe PMC article naming the gene
anywhere whose title or abstract names the trait. Accepted Novel|Hypothesized.

Selection signals are absent: the GWAS Catalog does not carry them, and no
comparable independent source was found. Sampling is seeded per trait, so the table can
be rebuilt; the catalog and both indexes move, so a rebuild months later may
not return the same pairs.

  python3 scripts/139_build_bench_bodyonly.py
"""
from __future__ import annotations

import csv
import os
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))
from pesto.services.pubmed_service import make_api_request_with_retry, NCBI_API_KEY  # noqa: E402
from pesto.config import BASE_URL_NCBI  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
HGNC = ROOT / "data" / "hgnc_complete_set.txt"
# The full associations table, from
# https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip
# (release of 15 September 2026). The REST API took 14 s a page and stopped
# answering deep pages of the large traits under load, so the catalog is read
# locally and only PubMed and Europe PMC are asked over the network.
CATALOG = ROOT / "data" / "gwas_catalog" / "gwas-catalog-download-associations-alt-full.tsv"
OUT = ROOT / "results" / "bench_bodyonly_labels.tsv"
BENCH = ROOT / "results" / "bench_bodyonly.tsv"
# Tables whose genes are excluded from the draw. The published draw excluded
# the genes of an earlier, unpublished test set, which is not distributed, so a
# rebuild from here can return other pairs.
EXISTING = []

# (catalog trait, pesto phenotype, type, PubMed title/abstract words)
TRAITS = [
    ("bipolar disorder", "bipolar disorder", "disease", ["bipolar"]),
    ("schizophrenia", "schizophrenia", "disease", ["schizophrenia", "schizophrenic"]),
    ("type 2 diabetes mellitus", "type 2 diabetes", "disease",
     ["type 2 diabetes", "T2D", "diabetes"]),
    ("asthma", "asthma", "disease", ["asthma"]),
    ("atrial fibrillation", "atrial fibrillation", "disease", ["atrial fibrillation"]),
    ("body height", "adult height", "trait", ["height", "stature"]),
    ("body mass index", "body mass index", "trait",
     ["body mass index", "BMI", "obesity", "adiposity"]),
    ("systolic blood pressure", "blood pressure", "trait",
     ["blood pressure", "hypertension"]),
]
POSITIVES_PER_TRAIT = 2
NEGATIVE_TRAITS = 6
SEED = 139
# Genes checked at once within a trait; eight traits run side by side, so
# about sixteen requests are in flight at the peak. Thirty-two made Europe
# PMC answer empty often enough to lose five traits of eight.
BATCH = 2


class Unanswered(Exception):
    """A service failed on every try. Under load Europe PMC and the catalog
    answer with an empty body, which reads exactly like "no hits", so a failure
    must not be allowed to pass for an answer."""


def get_json(url, params, tries=5, need=None):
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=90)
            if r.ok:
                d = r.json()
                if need is None or need in d:
                    return d
        except (requests.RequestException, ValueError):
            pass
        time.sleep(2 * (i + 1))
    raise Unanswered(url)


def hgnc():
    """Protein-coding symbols and their aliases."""
    names = {}
    with HGNC.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["locus_group"] != "protein-coding gene":
                continue
            alias = [a for a in (r["alias_symbol"] + "|" + r["prev_symbol"]).split("|")
                     if len(a) > 3 and a.isascii()]
            names[r["symbol"]] = alias[:4]
    return names


def already_used():
    used = set()
    for p in EXISTING:
        with p.open(encoding="utf-8") as fh:
            used |= {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}
    return used


def pubmed_count(gene, aliases, words):
    g = " OR ".join(f'"{n}"[tiab]' for n in [gene] + aliases)
    w = " OR ".join(f'"{x}"[tiab]' for x in words)
    params = {"db": "pubmed", "term": f"({g}) AND ({w})", "retmode": "json", "retmax": 0}
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", params)
    try:
        return int(r.json()["esearchresult"]["count"])
    except Exception:
        raise Unanswered("pubmed")


def epmc_count(query):
    return get_json(EPMC, {"query": query, "format": "json", "pageSize": 1},
                    need="hitCount")["hitCount"]


def named_in_body(gene, pmid):
    return epmc_count(f'EXT_ID:{pmid} AND SRC:MED AND '
                      f'(BODY:"{gene}" OR TABLE:"{gene}" OR SUPPL:"{gene}")') == 1


def load_catalog():
    """trait -> gene -> PMIDs, for associations at p <= 5e-8 lying inside one gene."""
    wanted = {t[0] for t in TRAITS}
    out = {t: {} for t in wanted}
    csv.field_size_limit(10 ** 7)
    with CATALOG.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            traits = wanted.intersection(t.strip() for t in r["MAPPED_TRAIT"].split(","))
            if not traits:
                continue
            gene = r["MAPPED_GENE"].strip()
            if not gene or r["INTERGENIC"] != "0" or any(c in gene for c in ",;- "):
                continue
            try:
                if float(r["P-VALUE"]) > 5e-8:
                    continue
            except ValueError:
                continue
            for t in traits:
                out[t].setdefault(gene, set()).add(r["PUBMEDID"].strip())
    return out


def candidates(trait, rng):
    order = sorted(CAT[trait])
    rng.shuffle(order)
    return order


def catalog_pmids(gene, trait):
    return CAT[trait].get(gene, set())


def check_positive(gene, trait, names, words):
    """Catalog articles naming the gene in their body, if the abstract route is blind."""
    try:
        if pubmed_count(gene, names[gene], words) != 0:
            return None
        pmids = sorted(catalog_pmids(gene, trait))
        body = []
        for p in pmids:
            if named_in_body(gene, p):
                body.append(p)
                if len(body) == 3:
                    break
    except Unanswered as e:
        print(f"! {gene:10} {trait}: no answer from {e}", flush=True)
        return None
    return (pmids, body) if body else None


def positives(trait, pheno, kind, words, names, used):
    # One generator per trait, so the pairs do not depend on which trait a
    # thread finished first. Candidates are checked in batches, and the first
    # POSITIVES_PER_TRAIT that pass in candidate order are kept.
    rng = random.Random(f"{SEED}:{trait}")
    order = [g for g in candidates(trait, rng) if g in names and g not in used]
    rows = []
    with ThreadPoolExecutor(BATCH) as ex:
        for i in range(0, len(order), BATCH):
            chunk = order[i:i + BATCH]
            for gene, hit in zip(chunk, ex.map(
                    lambda g: check_positive(g, trait, names, words), chunk)):
                if hit is None or len(rows) >= POSITIVES_PER_TRAIT:
                    continue
                pmids, body = hit
                rows.append({"gene": gene, "phenotype": pheno, "type": kind,
                             "accept": "Existing|Established",
                             "key_pmids": ";".join(body[:3]),
                             "status": "gwas_body_only",
                             "note": f"{len(pmids)} catalog GWAS article(s), "
                                     f"{len(body)} naming the gene in body/table/supplement",
                             "split": "bodyonly"})
                print(f"+ {gene:10} {pheno:22} {';'.join(body[:3])}", flush=True)
            if len(rows) >= POSITIVES_PER_TRAIT:
                break
    return rows


def check_negative(gene, trait, names, words):
    try:
        if pubmed_count(gene, names[gene], words) != 0:
            return False
        if catalog_pmids(gene, trait):
            return False
        where = " OR ".join(f'{f}:"{x}"' for f in ("TITLE", "ABSTRACT") for x in words)
        return epmc_count(f'"{gene}" AND ({where})') == 0
    except Unanswered as e:
        print(f"! {gene:10} {trait}: no answer from {e}", flush=True)
        return False


def negative(trait, pheno, kind, words, names, used):
    rng = random.Random(f"{SEED}:negative:{trait}")
    pool = [g for g in sorted(names) if g not in used]
    rng.shuffle(pool)
    with ThreadPoolExecutor(BATCH) as ex:
        for i in range(0, len(pool), BATCH):
            chunk = pool[i:i + BATCH]
            for gene, ok in zip(chunk, ex.map(
                    lambda g: check_negative(g, trait, names, words), chunk)):
                if ok:
                    print(f"- {gene:10} {pheno}", flush=True)
                    return [{"gene": gene, "phenotype": pheno, "type": kind,
                             "accept": "Novel|Hypothesized", "key_pmids": "",
                             "status": "negative_control",
                             "note": "no catalog association, no PubMed or Europe PMC pairing",
                             "split": "bodyonly"}]
    return []


CAT = {}


def main():
    CAT.update(load_catalog())
    for t in CAT:
        print(f"  {t}: {len(CAT[t])} genes", flush=True)
    names = hgnc()
    used = already_used()
    with ThreadPoolExecutor(len(TRAITS)) as ex:
        per_trait = list(ex.map(lambda t: positives(*t, names, used), TRAITS))
    rows, seen = [], set(used)
    for part in per_trait:
        for r in part:
            if r["gene"] not in seen:
                seen.add(r["gene"])
                rows.append(r)
    neg_traits = random.Random(SEED).sample(TRAITS, NEGATIVE_TRAITS)
    with ThreadPoolExecutor(len(neg_traits)) as ex:
        for part in ex.map(lambda t: negative(*t, names, seen), neg_traits):
            for r in part:
                if r["gene"] not in seen:
                    seen.add(r["gene"])
                    rows.append(r)

    cols = ["gene", "phenotype", "type", "accept", "key_pmids", "status", "note", "split"]
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with BENCH.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "phenotype"])
        w.writerows([r["gene"], r["phenotype"]] for r in rows)
    print(f"{len(rows)} pairs -> {OUT.relative_to(ROOT)}", file=sys.stderr)


if __name__ == "__main__":
    main()
