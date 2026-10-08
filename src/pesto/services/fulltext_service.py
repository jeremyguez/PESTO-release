"""Europe PMC: searching article bodies, and reading the part that names the gene.

PubMed indexes titles, abstracts and MeSH. A gene that appears only in a
results table, a supplementary file or one sentence of a discussion is
invisible to it, and that is not a rare corner: a genome-wide association
study often names its genes only in a results table or a supplement. A PubMed
search cannot reach those articles, however worded.

Two functions, and they are meant to be used together:

  search_fulltext   the gene in BODY/TABLE/SUPPL, the phenotype in the title
                    or abstract, which is the pairing PubMed cannot express
  gene_excerpts     the sentences around the gene in the full text, so that
                    an article whose abstract never mentions it still arrives
                    at the model with the reason it was retrieved

The excerpt matters as much as the search. Bands and the reading prompt see
title and abstract only, so a body-only article with an empty abstract would
be banded "no bearing on the question" and dropped before anything read it.

Everything here is free and needs no key. Failures return empty rather than
raising: a pipeline that already has a PubMed corpus should degrade to it
rather than stop.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import requests

SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

TIMEOUT = 45
# Europe PMC answers a malformed or throttled query with a bare {"version"},
# no hitCount and no error. Seen on roughly one call in thirty, and it looks
# exactly like a legitimate empty result, so every query is retried once.
RETRIES = 2


def _get(url, params=None, tries=RETRIES):
    for _ in range(tries):
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 200:
                return r
        except requests.RequestException:
            pass
    return None


def _quoted(terms):
    out, seen = [], set()
    for t in terms:
        t = (t or "").strip().replace('"', "")
        if len(t) > 2 and t.lower() not in seen:
            seen.add(t.lower())
            out.append(f'"{t}"')
    return out


def fulltext_query(gene, gene_aliases, phenotype, synonyms, anchors=()):
    """The gene in the body, the phenotype in the title or abstract.

    Aliases are capped: a gene with thirty published names makes a query Europe
    PMC rejects, and the long tail of those names is where the false hits are.
    Anchors are short terms any one of which places the phenotype, alongside
    the exact phrases rather than instead of them.
    """
    genes = _quoted([gene] + list(gene_aliases or [])[:4])
    pheno = _quoted([phenotype] + list(synonyms or [])[:8] + list(anchors or [])[:4])
    if not genes or not pheno:
        return ""
    body = " OR ".join(f"{f}:{g}" for f in ("BODY", "TABLE", "SUPPL") for g in genes)
    where = " OR ".join(f"{f}:{p}" for f in ("TITLE", "ABSTRACT") for p in pheno)
    return f"({body}) AND ({where})"


def search_fulltext(gene, gene_aliases, phenotype, synonyms, limit=25, anchors=()):
    """Articles whose body names the gene and whose abstract names the disease.

    Returns dicts in the shape the PubMed service returns, so a Corpus can be
    built from them without a special case, plus `pmcid`, which gene_excerpts
    needs and PubMed does not carry.
    """
    query = fulltext_query(gene, gene_aliases, phenotype, synonyms, anchors)
    if not query:
        return []
    r = _get(SEARCH, {"query": query, "format": "json", "pageSize": limit,
                      "resultType": "lite"})
    if r is None:
        return []
    data = r.json()
    if "hitCount" not in data:
        r = _get(SEARCH, {"query": query, "format": "json", "pageSize": limit,
                          "resultType": "lite"}, tries=1)
        data = r.json() if r is not None else {}
    out = []
    for a in (data.get("resultList") or {}).get("result") or []:
        pmid = (a.get("pmid") or "").strip()
        if not pmid:
            continue
        out.append({"pmid": pmid, "pmcid": a.get("pmcid") or "",
                    "title": a.get("title") or "",
                    "abstract": "", "pubdate": str(a.get("pubYear") or ""),
                    "query_tier": "fulltext"})
    return out


def pmcids_for(pmids):
    """PMC identifiers for PubMed ids, in one query per fifty.

    The Corpus type carries no pmcid, and adding a field to it would touch
    every run ever stored, so the id is looked up again when an excerpt is
    wanted. Articles with no open-access copy simply do not come back.
    """
    out = {}
    pmids = [p for p in pmids if p]
    for i in range(0, len(pmids), 50):
        chunk = pmids[i:i + 50]
        # The parentheses are load-bearing: Europe PMC binds AND tighter than
        # OR, so without them only the last id in the chunk is constrained to
        # SRC:MED and the rest come back empty.
        query = "(" + " OR ".join(f"EXT_ID:{p}" for p in chunk) + ") AND SRC:MED"
        r = _get(SEARCH, {"query": query, "format": "json",
                          "pageSize": len(chunk), "resultType": "lite"})
        if r is None:
            continue
        for a in ((r.json().get("resultList") or {}).get("result") or []):
            if a.get("pmcid"):
                out[(a.get("pmid") or "").strip()] = a["pmcid"]
    return out


def _plain(xml_text):
    """Body text with tables flattened, tags dropped, whitespace collapsed."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml_text))
    parts = []
    for el in root.iter():
        if el.tag in ("td", "th"):
            parts.append((el.text or "").strip() + " |")
        elif el.text:
            parts.append(el.text)
        if el.tail:
            parts.append(el.tail)
    return re.sub(r"\s+", " ", " ".join(parts))


def fetch_fulltext(pmcid):
    """Open-access full text, Europe PMC first and NCBI as the fallback.

    Author manuscripts are the reason for the second call: Europe PMC answers
    403 for some of them (PMC5600716 among others) while efetch serves them.
    """
    if not pmcid:
        return ""
    r = _get(FULLTEXT.format(pmcid=pmcid), tries=1)
    if r is not None and r.text.lstrip().startswith("<"):
        return _plain(r.text)
    r = _get(EFETCH, {"db": "pmc", "id": pmcid.replace("PMC", ""),
                      "retmode": "xml"}, tries=1)
    return _plain(r.text) if r is not None else ""


def gene_excerpts(pmcid, gene, gene_aliases=(), window=420, most=3):
    """Up to `most` passages around the gene, longest-name match winning.

    Aliases are searched after the symbol and only if it is absent, because a
    two-letter alias matches everything: on a paper that writes both, the
    symbol is what the authors meant.
    """
    text = fetch_fulltext(pmcid)
    if not text:
        return []
    names = [gene] + [a for a in (gene_aliases or []) if len(a) > 3]
    spans = []
    for name in names:
        spans = [m.start() for m in
                 re.finditer(rf"\b{re.escape(name)}\b", text, re.I)]
        if spans:
            break
    out = []
    for pos in spans[:most]:
        start, end = max(0, pos - window // 2), pos + window // 2
        piece = text[start:end].strip()
        if start:
            piece = "..." + piece
        if end < len(text):
            piece += "..."
        out.append(piece)
    return out
