"""The Open Targets branch: a gene's associated traits, and whether the
queried phenotype is among them.

The gene symbol is resolved to an Ensembl id, every associated disease or trait
is fetched from the Open Targets Platform GraphQL API, those scoring above 0.05
are handed to ot_matcher, and the trait it names (if any) decides the verdict
in scoring.open_targets_verdict. The matcher's record is saved beside the run,
keyed by its settings, so asking again reuses it without a model call.
"""
import json
import os
import re

import requests

from . import ot_matcher
from .config import BASE_URL_NCBI
from .services.pubmed_service import NCBI_API_KEY, make_api_request_with_retry

OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def _run_open_targets_query(query, variables=None):
    response = requests.post(
        OPEN_TARGETS_URL,
        json={"query": query, "variables": variables or {}},
        timeout=20
    )
    response.raise_for_status()
    payload = response.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    return payload.get("data")


def _get_ensembl_id_from_symbol(gene_symbol):
    query = """
    query searchGene($q: String!) {
      search(queryString: $q, entityNames: ["target"]) {
        hits {
          id
          name
          __typename
        }
      }
    }
    """
    data = _run_open_targets_query(query, {"q": gene_symbol})
    hits = (data or {}).get("search", {}).get("hits") or []
    if not hits:
        print(f"ERROR: Open Targets search returned no hits for '{gene_symbol}'. Raw data: {data}")
        return None
    # Pick first hit; prefer exact symbol match if present
    gene_upper = (gene_symbol or "").upper()
    exact = next((h for h in hits if (h.get("name") or "").upper() == gene_upper), None)
    chosen = exact or hits[0]
    ensembl_id = chosen.get("id")
    if not ensembl_id:
        print(f"ERROR: No Ensembl ID found in hit for '{gene_symbol}': {chosen}")
        return None
    return ensembl_id


def _extract_pmid(lit_entry):
    """Normalize a literature entry to a PMID string."""
    if isinstance(lit_entry, dict):
        return lit_entry.get("pmid") or lit_entry.get("pmId") or lit_entry.get("id")
    if lit_entry is None:
        return None
    return str(lit_entry)


def _enrich_with_pubmed_metadata(articles):
    """Fetch titles/pubdates from PubMed for given PMIDs and enrich in place."""
    if not articles:
        return articles
    pmids = [a["pmid"] for a in articles if a.get("pmid")]
    if not pmids:
        return articles

    esummary_url = f"{BASE_URL_NCBI}esummary.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json"
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    resp = make_api_request_with_retry(esummary_url, params)
    if not resp:
        return articles

    try:
        summary = resp.json().get("result", {})
    except Exception as exc:
        print(f"ERROR: Failed to parse PubMed esummary response: {exc}")
        return articles

    lookup = {a["pmid"]: a for a in articles}
    for pmid in pmids:
        info = summary.get(pmid) or {}
        art = lookup.get(pmid)
        if not art:
            continue
        title = info.get("title")
        pubdate = info.get("pubdate")
        if title:
            art["title"] = title
        if pubdate:
            art["year"] = pubdate
    return articles


def _fetch_evidence_for_target_disease(ensembl_id, disease_id, size=25):
    """Fetch literature evidence from Open Targets, then enrich with PubMed titles."""
    query = """
    query targetEvidence($ensemblId: String!, $efoId: String!, $size: Int!) {
      target(ensemblId: $ensemblId) {
        evidences(efoIds: [$efoId], size: $size) {
          rows {
            literature
            publicationYear
            datasourceId
          }
        }
      }
    }
    """
    variables = {"ensemblId": ensembl_id, "efoId": disease_id, "size": size}
    data = _run_open_targets_query(query, variables)
    rows = (
        (data or {})
        .get("target", {})
        .get("evidences", {})
        .get("rows", [])
        or []
    )
    if not rows:
        return []

    articles = []
    seen_pmids = set()
    for row in rows:
        lit_list = row.get("literature") or []
        year = row.get("publicationYear")
        source = row.get("datasourceId")
        for lit in lit_list:
            pmid = _extract_pmid(lit)
            if not pmid or pmid in seen_pmids:
                continue
            seen_pmids.add(pmid)
            articles.append({
                "pmid": str(pmid),
                "title": "",  # will be filled from PubMed
                "abstract": "",
                "year": year,
                "source": source
            })

    return _enrich_with_pubmed_metadata(articles)


def _derive_open_targets_verdict(open_targets_traits):
    """The Open Targets verdict. The rule lives in scoring.py."""
    from .scoring import open_targets_verdict
    return open_targets_verdict(open_targets_traits)


def _fetch_traits_for_ensembl_id(ensembl_id, page_size=200):
    # Map disease_id -> {"name": str, "score": float}
    traits = {}
    page_index = 0
    query = """
    query targetDiseaseAssociations($ensembl_id: String!, $size: Int!, $index: Int!) {
      target(ensemblId: $ensembl_id) {
        associatedDiseases(page: { size: $size, index: $index }, includeMeasurements: true) {
          count
          rows {
            score
            disease {
              id
              name
            }
          }
        }
      }
    }
    """
    while True:
        variables = {
            "ensembl_id": ensembl_id,
            "size": page_size,
            "index": page_index
        }
        data = _run_open_targets_query(query, variables)
        target = (data or {}).get("target")
        if not target:
            break
        assoc = target.get("associatedDiseases") or {}
        rows = assoc.get("rows") or []
        total_count = assoc.get("count") or 0
        for row in rows:
            disease = row.get("disease") or {}
            disease_id = disease.get("id")
            disease_name = disease.get("name")
            score = row.get("score")
            if disease_id and disease_name:
                existing = traits.get(disease_id)
                # Keep max score if multiple rows exist
                if not existing or (score is not None and score > existing.get("score", float("-inf"))):
                    traits[disease_id] = {"name": disease_name, "score": score}
        if (page_index + 1) * page_size >= total_count or not rows:
            break
        page_index += 1
    return traits


def _sanitize_gene_for_filename(gene_symbol):
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene_symbol or "")
    return sanitized or "gene"


def _ot_match_sidecar_path(runs_dir, gene_symbol, phenotype):
    return os.path.join(
        runs_dir,
        "open_targets_match_{}_{}.json".format(
            _sanitize_gene_for_filename(gene_symbol),
            _sanitize_gene_for_filename(phenotype)))


def _traits_txt_path(runs_dir, gene_symbol, phenotype):
    return os.path.join(
        runs_dir,
        "open_targets_traits_{}_{}.txt".format(
            _sanitize_gene_for_filename(gene_symbol),
            _sanitize_gene_for_filename(phenotype)))


def _write_traits_txt(path, traits):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for t in traits:
            score = t.get("score")
            score_str = "" if score is None else str(score)
            f.write("{}\t{}\t{}\t{}\n".format(
                t.get("disease_id", ""), t.get("name", ""), score_str,
                bool(t.get("matched"))))


def _apply_match_to_traits(traits, match_detail):
    matched_id = match_detail.get("matched_id")
    out = []
    for t in traits:
        row = dict(t)
        row["matched"] = bool(matched_id) and row.get("disease_id") == matched_id
        out.append(row)
    return out


def _result_from_record(record, gene_symbol, phenotype, runs_dir, gate=None):
    """Rebuild the fetch payload from a saved agent record. No model call."""
    derived = ot_matcher.replay_tagged(record, apply_gate=gate)
    traits = _apply_match_to_traits(record.get("traits") or [], derived)
    txt = _traits_txt_path(runs_dir, gene_symbol, phenotype)
    _write_traits_txt(txt, traits)
    match = {k: v for k, v in derived.items()}
    match.update({
        "shortlist": record.get("shortlist") or [],
        "gate_text": record.get("gate_text") or "",
        "gate_dropped_ids": record.get("gate_dropped_ids") or [],
        "tag_text": record.get("tag_text") or "",
        "graded": record.get("graded") or [],
        "tokens": record.get("tokens") or {},
        "cached": True,
    })
    return {
        "ensembl_id": record.get("ensembl_id"),
        "traits_count": len(traits),
        "file_path": txt,
        "traits": traits,
        "articles": record.get("articles") or [],
        "ot_synonym_selection_usage": None,
        "match": match,
        "cached": True,
    }


def _load_match_sidecar(path, matcher=None, encoder=None, top_k=None,
                        reasons=None):
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            record = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None
    from . import config
    want = config.ot_agent_version(matcher, encoder, top_k, reasons)
    if record.get("agent_version") != want:
        return None
    if not record.get("shortlist") or not record.get("graded"):
        return None
    return record


def _write_match_sidecar(path, record):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, default=str)
        fh.write("\n")


def fetch_and_save_open_targets_traits(gene_symbol, runs_dir, phenotype,
                                       synonyms, model="claude-haiku-4-5",
                                       temperature=0.0, tracker=None,
                                       matcher=None, encoder=None, top_k=None,
                                       gate=None, reasons=None):
    sidecar = _ot_match_sidecar_path(runs_dir, gene_symbol, phenotype)
    cached = _load_match_sidecar(sidecar, matcher, encoder, top_k, reasons)
    if cached is not None:
        print(f"INFO: Open Targets - reused {sidecar}")
        return _result_from_record(cached, gene_symbol, phenotype,
                                   runs_dir, gate=gate)
    try:
        ensembl_id = _get_ensembl_id_from_symbol(gene_symbol)
        if not ensembl_id:
            result = {
                "error": f"No Ensembl ID found for symbol '{gene_symbol}'",
                "stage": "search"
            }
            print(f"INFO: Open Targets - {result}")
            return result
        traits = _fetch_traits_for_ensembl_id(ensembl_id)
        os.makedirs(runs_dir, exist_ok=True)
        file_path = _traits_txt_path(runs_dir, gene_symbol, phenotype)
        # Preserve insertion order, filter score > 0.05
        filtered_items = []
        for disease_id, payload in traits.items():
            score = payload.get("score")
            if score is None or score <= 0.05:
                continue
            filtered_items.append((disease_id, payload))

        # Which trait, if any, records the reported association. The matcher
        # defaults to config.OT_MATCHER and names at most one trait, so the
        # verdict below rests on a trait that can be pointed at rather than on
        # whichever of a dozen loose string hits happened to score highest.
        # `synonyms` is no longer consulted: no matcher searches for a literal
        # string, which is what the synonym list was for.
        match_detail = ot_matcher.match(gene_symbol, phenotype,
                                        [{"disease_id": did,
                                          "name": payload.get("name", ""),
                                          "score": payload.get("score")}
                                         for did, payload in filtered_items],
                                        method=matcher, model=None,
                                        encoder=encoder, top_k=top_k,
                                        gate=gate, reasons=reasons)
        labels_map = {did: did == match_detail.get("matched_id")
                      for did, _ in filtered_items}
        ot_synonym_selection_usage = match_detail.get("usage")

        traits_payload = []
        with open(file_path, "w", encoding="utf-8") as f:
            for disease_id, payload in filtered_items:
                name = payload.get("name", "")
                score = payload.get("score")
                score_str = "" if score is None else str(score)
                is_match = labels_map.get(disease_id, False)
                f.write(f"{disease_id}\t{name}\t{score_str}\t{is_match}\n")
                traits_payload.append({
                    "disease_id": disease_id,
                    "name": name,
                    "score": score,
                    "matched": is_match
                })
        # Collect matched traits sorted by score desc
        matched_sorted = sorted(
            [t for t in traits_payload if t["matched"] and t["score"] is not None],
            key=lambda x: x["score"],
            reverse=True
        )
        ot_articles = []
        for t in matched_sorted:
            if len(ot_articles) >= 10:
                break
            # Fetch evidence for each matched trait (don't let this fail the whole process)
            try:
                articles = _fetch_evidence_for_target_disease(ensembl_id, t["disease_id"], size=25)
                for art in articles:
                    if len(ot_articles) >= 10:
                        break
                    ot_articles.append(art)
            except Exception as e:
                print(f"WARNING: Failed to fetch evidence for {t['disease_id']}: {e}")
        from . import config
        match_public = {k: v for k, v in match_detail.items() if k != "usage"}
        record = {
            "gene": gene_symbol,
            "phenotype": phenotype,
            "ensembl_id": ensembl_id,
            "agent_version": config.ot_agent_version(
                matcher, encoder, top_k, reasons),
            "ot_match_version": config.ot_match_version(
                matcher, encoder, top_k, gate, reasons),
            "traits": [{"disease_id": t["disease_id"],
                        "name": t["name"],
                        "score": t["score"]}
                       for t in traits_payload],
            "shortlist": match_detail.get("shortlist") or [],
            "gate_text": match_detail.get("gate_text") or "",
            "gate_dropped_ids": match_detail.get("gate_dropped_ids") or [],
            "tag_text": match_detail.get("tag_text") or "",
            "graded": match_detail.get("graded") or [],
            "of": match_detail.get("of"),
            "tokens": match_detail.get("tokens") or {},
            "articles": ot_articles,
        }
        _write_match_sidecar(sidecar, record)
        result = {
            "ensembl_id": ensembl_id,
            "traits_count": len(filtered_items),
            "file_path": file_path,
            "traits": traits_payload,
            "articles": ot_articles,
            "ot_synonym_selection_usage": ot_synonym_selection_usage,
            "match": match_public,
        }
        print(f"INFO: Open Targets - success {result}")
        return result
    except Exception as e:
        error_msg = f"Open Targets fetch failed for {gene_symbol}: {e}"
        result = {
            "error": error_msg,
            "stage": "fetch"
        }
        print(f"ERROR: {result}")
        return result
