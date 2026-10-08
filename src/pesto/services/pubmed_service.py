"""PubMed API client for searching and fetching articles."""
import random
import re
import threading
import time
import requests
from xml.etree import ElementTree
from datetime import datetime
from ..config import BASE_URL_NCBI
from ..credentials import ncbi_api_key

# Optional. NCBI serves unkeyed callers too, at three requests a second instead
# of ten, and the rate gate below reads this to pick its pace. Nothing is
# shipped here: a key identifies whoever registered it, so it belongs to the
# person running the pipeline, not to the code.
NCBI_API_KEY = ncbi_api_key()
PUBMED_CALL_COOLDOWN = 0.05
VERBOSE_MODE = False  # Global flag to control verbose output

# Optional publication-date ceiling for ALL PubMed esearch queries. When set
# (e.g. "2019/12/31"), every search is restricted to articles published on or
# before this date via datetype=pdat. Default None = no restriction (historical
# behavior unchanged). Used to reconstruct the literature as it stood at a past
# date (e.g. a "2019" run to test prospective discovery potential).
MAX_PUBDATE = None
MIN_PUBDATE = "1000/01/01"


def _apply_pubdate_cap(esearch_params):
    """Restrict an esearch params dict to MAX_PUBDATE if one is configured.

    Only esearch needs this; esummary/efetch fetch by explicit PMID so they
    inherit the date restriction for free. Mutates and returns the dict.
    """
    if MAX_PUBDATE:
        esearch_params["datetype"] = "pdat"
        esearch_params["mindate"] = MIN_PUBDATE
        esearch_params["maxdate"] = MAX_PUBDATE
    return esearch_params


def _filter_articles_by_pubdate_cap(articles):
    """Deterministically drop any article whose displayed publication year is
    after the configured cap.

    The server-side datetype=pdat restriction occasionally lets through records
    whose electronic pub date is <= cap but whose journal-issue date (the value
    shown in `pubdate`) is later. For a clean, defensible "as of <year>" run we
    also enforce the ceiling on the visible year here.
    """
    if not MAX_PUBDATE:
        return articles
    cap_year = int(str(MAX_PUBDATE)[:4])
    kept = []
    for a in articles:
        m = re.search(r"(19|20)\d{2}", str(a.get("pubdate", "")))
        # Drop only when we can parse a year AND it exceeds the cap. Unparseable
        # dates are kept (server-side pdat already restricted them).
        if m and int(m.group(0)) > cap_year:
            continue
        kept.append(a)
    return kept

# Gènes avec noms ambigus qui nécessitent une recherche PubMed spécifique
# Pour ces gènes, on utilise "X gene"[Title/Abstract] au lieu de "X"[Title/Abstract]
# pour éviter le bruit (mots anglais courants, abréviations médicales, noms de maladies)
AMBIGUOUS_GENE_NAMES = {
    # Abréviations médicales/scientifiques
    "APC", "BAD", "BAX", "BID", "CA2", "CAMP", "EGF", "EGFR", "FOS", "GC", "HR",
    "JUN", "MTOR", "MYC", "NGF", "PIP", "SRC", "TNF", "VHL",
    # Noms de maladies/syndromes
    "CAD", "DMD", "FAP", "SCD",
    # Mots anglais courants en science
    "APP", "ARC", "CAT", "IMPACT", "KIN", "KIT", "MAG", "MAX", "MET",
    "NODAL", "RAN", "REST", "SET", "SON", "TANK", "TUB", "VIM",
    # Gènes de 2 caractères (trop ambigus)
    "AR", "C2", "C3", "C5", "C6", "C7", "C9", "CP", "CS", "F2", "F3", "F5",
    "F7", "F8", "F9", "FH", "GK", "HP", "IK", "KL", "KY", "MB", "PC", "SI",
    "TF", "TG", "TH", "XG", "XK",
}


def get_pubmed_gene_term(gene_name: str) -> str:
    """
    Returns the appropriate PubMed search term for a gene.
    
    For ambiguous genes (common English words, medical abbreviations, disease names),
    uses "GENE gene"[Title/Abstract] to reduce noise.
    For other genes, uses "GENE"[Title/Abstract].
    
    Args:
        gene_name: The gene symbol (e.g., "SET", "BRCA1")
        
    Returns:
        The PubMed search term (e.g., '"SET gene"[Title/Abstract]')
    """
    if gene_name.upper() in AMBIGUOUS_GENE_NAMES:
        return f'"{gene_name} gene"[Title/Abstract]'
    else:
        return f'"{gene_name}"[Title/Abstract]'


# NCBI allows ten requests per second with a key and three without, counted
# across the whole client rather than per connection. Spacing requests inside
# each caller cannot enforce that: it is either too slow when one thread runs,
# or too fast when twelve do. One shared gate does both.
NCBI_MAX_PER_SECOND = 9.0 if NCBI_API_KEY else 2.5
_rate_gate = threading.Lock()
_next_slot = [0.0]


def _wait_for_slot():
    with _rate_gate:
        now = time.monotonic()
        slot = max(now, _next_slot[0])
        _next_slot[0] = slot + 1.0 / NCBI_MAX_PER_SECOND
    if slot > now:
        time.sleep(slot - now)


def make_api_request_with_retry(url, params, timeout=30):
    """Makes an API request with an extended retry mechanism."""
    delays = [0.1] * 10 + [0.5] * 10 + [1] * 10 + [2, 5, 10, 20, 30, 60, 120]
    max_attempts = len(delays) + 1
    for i in range(max_attempts):
        try:
            _wait_for_slot()
            if PUBMED_CALL_COOLDOWN > 0:
                time.sleep(PUBMED_CALL_COOLDOWN)
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response
        except (requests.exceptions.RequestException, requests.exceptions.HTTPError) as e:
            if VERBOSE_MODE:
                print(f"WARNING: API request failed (attempt {i+1}/{max_attempts}): {e}")
            if i < len(delays):
                time.sleep(delays[i])
            else:
                # Toujours afficher les erreurs critiques
                if not VERBOSE_MODE:
                    print("ERROR: All API request attempts failed.")
                else:
                    print("ERROR: All API request attempts failed.")
                return None
    return None


def fetch_abstracts(pmids):
    """
    Fetch abstracts for a list of PMIDs using efetch.
    Returns a dictionary mapping PMID -> abstract text.
    """
    if not pmids:
        return {}
    
    efetch_url = f"{BASE_URL_NCBI}efetch.fcgi"
    
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "rettype": "abstract"
    }
    
    if NCBI_API_KEY:
        params['api_key'] = NCBI_API_KEY
    
    response = make_api_request_with_retry(efetch_url, params)
    if not response:
        return {}
    
    abstracts = {}
    
    try:
        root = ElementTree.fromstring(response.content)
        
        for article in root.findall('.//PubmedArticle'):
            pmid_elem = article.find('.//PMID')
            # Every <AbstractText>, not the first. Structured abstracts
            # (BACKGROUND / METHODS / RESULTS / CONCLUSIONS) carry one element
            # per section, and the results are rarely in the first.
            sections = []
            for el in article.findall('.//Abstract/AbstractText'):
                text = "".join(el.itertext()).strip()
                if not text:
                    continue
                label = (el.get('Label') or "").strip()
                sections.append(f"{label}: {text}" if label else text)
            
            if pmid_elem is not None:
                abstracts[pmid_elem.text] = " ".join(sections)
        
        if VERBOSE_MODE:
            print(f"INFO: Fetched abstracts for {len(abstracts)}/{len(pmids)} articles.")
        
    except Exception as e:
        print(f"ERROR: Could not parse abstract XML: {e}")
    
    return abstracts


def search_pubmed(gene_name, keywords, max_results):
    """
    Searches PubMed using a single combined query (always one-search mode).
    Deduplicates results and tracks how each article was found.
    """
    if VERBOSE_MODE:
        print("INFO: Performing combined single search.")
    
    # Use appropriate search term based on gene name ambiguity
    gene_term = get_pubmed_gene_term(gene_name)
    
    if keywords:
        keyword_clause = " OR ".join([f'"{kw}"' for kw in keywords])
        search_term = f'({gene_term}) AND ({keyword_clause})'
    else:
        search_term = gene_term
    
    search_queries = [{"term": search_term, "tag": "one_search", "retmax": max_results}]

    all_results = {}

    for query in search_queries:
        term, tag, retmax = query["term"], query["tag"], query["retmax"]

        esearch_params = {
            "db": "pubmed", "term": term, "retmax": retmax,
            "sort": "relevance", "usehistory": "y", "retmode": "json"
        }
        _apply_pubdate_cap(esearch_params)
        if NCBI_API_KEY:
            esearch_params['api_key'] = NCBI_API_KEY

        search_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", esearch_params)
        if not search_response:
            continue

        try:
            id_list = search_response.json().get("esearchresult", {}).get("idlist", [])
            if not id_list:
                continue

            esummary_params = {"db": "pubmed", "id": ",".join(id_list), "retmode": "json"}
            if NCBI_API_KEY:
                esummary_params['api_key'] = NCBI_API_KEY

            summary_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", esummary_params)
            if not summary_response:
                continue
            
            summary_data = summary_response.json()
            
            for pmid in id_list:
                if pmid not in all_results:
                    article_info = summary_data.get("result", {}).get(pmid, {})
                    title = article_info.get("title", "No title found")
                    pubdate = article_info.get("pubdate", "No date")
                    all_results[pmid] = {
                        "title": title, 
                        "pmid": pmid, 
                        "pubdate": pubdate,
                        "found_by": [tag],
                        "abstract": None
                    }
                else:
                    if tag not in all_results[pmid]["found_by"]:
                        all_results[pmid]["found_by"].append(tag)
        
        except Exception as e:
            print(f"ERROR: Could not parse PubMed response for term '{term}': {e}")
            continue

    final_results = _filter_articles_by_pubdate_cap(list(all_results.values()))
    if VERBOSE_MODE:
        print(f"INFO: Found a total of {len(final_results)} unique articles for {gene_name}.")
    
    return final_results


def search_pubmed_query(search_term, max_results=50):
    """
    Executes a single PubMed query string and returns summary records.
    """
    if VERBOSE_MODE:
        print(f"INFO: Custom PubMed query: {search_term}")

    esearch_params = {
        "db": "pubmed",
        "term": search_term,
        "retmax": max_results,
        "sort": "relevance",
        "retmode": "json"
    }
    _apply_pubdate_cap(esearch_params)
    if NCBI_API_KEY:
        esearch_params["api_key"] = NCBI_API_KEY

    search_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", esearch_params)
    if not search_response:
        return []

    try:
        id_list = search_response.json().get("esearchresult", {}).get("idlist", [])
    except Exception as exc:
        print(f"ERROR: Failed to parse custom query esearch response: {exc}")
        return []

    if not id_list:
        return []

    esummary_params = {"db": "pubmed", "id": ",".join(id_list), "retmode": "json"}
    if NCBI_API_KEY:
        esummary_params["api_key"] = NCBI_API_KEY

    summary_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", esummary_params)
    if not summary_response:
        return []

    try:
        summary_data = summary_response.json()
    except Exception as exc:
        print(f"ERROR: Failed to parse custom query esummary response: {exc}")
        return []

    results = []
    for pmid in id_list:
        article_info = summary_data.get("result", {}).get(pmid, {})
        if not article_info:
            continue
        results.append(
            {
                "title": article_info.get("title", "No title found"),
                "pmid": pmid,
                "pubdate": article_info.get("pubdate", "No date"),
                "found_by": ["custom_query"],
                "abstract": None,
            }
        )

    return _filter_articles_by_pubdate_cap(results)


def search_lof_gof_pubmed(gene_name, max_results=10):
    """
    Search PubMed for LoF/GoF related articles for a gene.
    Returns a list of articles with pmid, title, and found_by tags.
    """
    esearch_url = f"{BASE_URL_NCBI}esearch.fcgi"
    esummary_url = f"{BASE_URL_NCBI}esummary.fcgi"
    
    # Use appropriate search term based on gene name ambiguity
    gene_term = get_pubmed_gene_term(gene_name)
    query = f"{gene_term} AND (gain-of-function OR loss-of-function)"
    
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "sort": "relevance"
    }
    _apply_pubdate_cap(params)
    
    if VERBOSE_MODE:
        print(f"INFO: Searching PubMed for LoF/GoF: {query}")
    
    search_response_raw = make_api_request_with_retry(esearch_url, params)
    if not search_response_raw:
        return []
    
    try:
        search_response = search_response_raw.json()
    except Exception as e:
        print(f"ERROR: Failed to parse search response JSON: {e}")
        return []
    
    pmids = search_response.get("esearchresult", {}).get("idlist", [])
    
    if not pmids:
        if VERBOSE_MODE:
            print(f"INFO: No LoF/GoF articles found for {gene_name}")
        return []
    
    if VERBOSE_MODE:
        print(f"INFO: Found {len(pmids)} LoF/GoF articles")
    
    summary_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json"
    }
    
    summary_response_raw = make_api_request_with_retry(esummary_url, summary_params)
    if not summary_response_raw:
        return []
    
    try:
        summary_response = summary_response_raw.json()
    except Exception as e:
        print(f"ERROR: Failed to parse summary response JSON: {e}")
        return []
    
    articles = []
    result_dict = summary_response.get("result", {})
    
    for pmid in pmids:
        if pmid in result_dict:
            article_data = result_dict[pmid]
            articles.append({
                "pmid": pmid,
                "title": article_data.get("title", "No title available"),
                "pubdate": article_data.get("pubdate", "No date"),
                "found_by": ["lof_gof_search"],
                "abstract": None
            })
    
    return _filter_articles_by_pubdate_cap(articles)


def search_pubmed_gene_phenotype(
    gene_name,
    phenotype,
    phenotype_terms,
    max_results_primary=30,
    max_results_secondary=20,
    delay_ms=0,
    gene_aliases=None,
    max_results_oldest=0,
    require_terms=None
):
    """
    Searches PubMed for articles mentioning:
    1) gene AND any phenotype term (primary search, max_results_primary)
    2) gene AND (phenotype OR disease OR syndrome OR association OR trait) (secondary search, max_results_secondary)
    3) the oldest hits of the primary term (max_results_oldest, off by default)

    `require_terms`, if given, is AND-ed onto the phenotype clause as a
    Title/Abstract disjunction, so a name search can ask for papers that
    also use the vocabulary of genetic evidence.

    The third pass exists because relevance ranking favours recent and heavily
    cited work, so for a gene settled twenty years ago it returns reviews and
    leaves the founding report out of the corpus. Reserving a few slots for the
    oldest hits puts the first description back in reach, and unlike relevance
    the tail of a date-sorted list is stable from one run to the next.

    Combines, deduplicates, then fetches summaries.
    Returns (articles, alias_hits) where alias_hits are aliases that contributed at least one unique PMID.
    """
    print(f"INFO: Searching PubMed for {gene_name} AND phenotypes: {phenotype_terms}")
    
    def _maybe_random_delay():
        if not delay_ms or delay_ms <= 0:
            return
        wait_ms = random.uniform(delay_ms, delay_ms * 10)
        time.sleep(wait_ms / 1000.0)

    def _run_esearch(term, retmax, sort="relevance", retstart=0):
        params = {
            "db": "pubmed",
            "term": term,
            "retmax": retmax,
            "retstart": retstart,
            "sort": sort,
            "usehistory": "y",
            "retmode": "json"
        }
        if NCBI_API_KEY:
            params['api_key'] = NCBI_API_KEY
        _maybe_random_delay()
        resp = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", params)
        _maybe_random_delay()
        return resp
    
    def _fetch_summaries(pmids):
        if not pmids:
            return None
        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json"
        }
        if NCBI_API_KEY:
            params['api_key'] = NCBI_API_KEY
        _maybe_random_delay()
        resp = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", params)
        _maybe_random_delay()
        return resp
    
    def _apply_random_delay():
        if not delay_ms or delay_ms <= 0:
            return
        wait_ms = random.uniform(delay_ms, delay_ms * 10)
        time.sleep(wait_ms / 1000.0)

    # Build clauses: one for the main gene, one OR-combined for aliases (without main)
    main_gene = (gene_name or "").strip()
    alias_tokens = []
    seen_tokens = set()
    for tok in list(gene_aliases or []):
        normalized = (tok or "").strip()
        if not normalized:
            continue
        upper = normalized.upper()
        if upper in seen_tokens or upper == main_gene.upper():
            continue
        seen_tokens.add(upper)
        alias_tokens.append(normalized)

    # Text Word also covers the substance names NLM assigns, so a paper is found
    # when the indexers tagged the gene but the authors only spelled the protein
    # out. Title/Abstract alone misses those; All Fields would add authors and
    # affiliations, which is noise for symbols that are also common words.
    gene_clause_main = f'"{main_gene}"[Text Word]' if main_gene else ""
    gene_clause_alias = " OR ".join([f'"{tok}"[Text Word]' for tok in alias_tokens])

    # Queries
    phenotype_query = " OR ".join([f'"{term}"[Title/Abstract]' for term in phenotype_terms])
    if require_terms:
        needed = " OR ".join([f'"{term}"[Title/Abstract]' for term in require_terms])
        phenotype_query = f"({phenotype_query}) AND ({needed})"
    secondary_query = '"phenotype"[Title/Abstract] OR "disease"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "association"[Title/Abstract] OR "trait"[Title/Abstract]'

    primary_pmids = []
    secondary_pmids = []
    seen_pmids = set()
    # The oldest-first pass needs the size of the same two result sets, to know
    # how far in to start reading. esearch reports it on every response, so the
    # figures are kept here rather than asked for a second time.
    primary_total = None
    secondary_total = None

    alias_hits = []

    def _collect(pmids, target_list, max_len):
        added = False
        for pmid in pmids:
            if pmid not in seen_pmids and len(target_list) < max_len:
                target_list.append(pmid)
                seen_pmids.add(pmid)
                added = True
        return added

    # Main gene searches
    try:
        if gene_clause_main:
            _apply_random_delay()
            if max_results_primary > 0:
                primary_term_main = f'(({gene_clause_main})) AND ({phenotype_query})'
                primary_resp_main = _run_esearch(primary_term_main, max_results_primary)
                if primary_resp_main:
                    result = primary_resp_main.json().get("esearchresult", {})
                    primary_total = int(result.get("count", 0))
                    pmids = result.get("idlist", [])
                    _collect(pmids, primary_pmids, max_results_primary)
            _apply_random_delay()
            if max_results_secondary > 0:
                secondary_term_main = f'(({gene_clause_main})) AND ({secondary_query})'
                secondary_resp_main = _run_esearch(secondary_term_main, max_results_secondary)
                if secondary_resp_main:
                    result = secondary_resp_main.json().get("esearchresult", {})
                    secondary_total = int(result.get("count", 0))
                    pmids = result.get("idlist", [])
                    _collect(pmids, secondary_pmids, max_results_secondary)
    except Exception as e:
        print(f"ERROR: Failed to parse PubMed search results for main gene: {e}")

    # Alias OR searches (excluding main gene)
    try:
        if gene_clause_alias:
            remaining_primary = max_results_primary - len(primary_pmids)
            remaining_secondary = max_results_secondary - len(secondary_pmids)

            if remaining_primary > 0:
                _apply_random_delay()
                primary_term_alias = f'(({gene_clause_alias})) AND ({phenotype_query})'
                primary_resp_alias = _run_esearch(primary_term_alias, remaining_primary)
                if primary_resp_alias:
                    pmids = primary_resp_alias.json().get("esearchresult", {}).get("idlist", [])
                    if _collect(pmids, primary_pmids, max_results_primary):
                        alias_hits = alias_tokens

            if remaining_secondary > 0:
                _apply_random_delay()
                secondary_term_alias = f'(({gene_clause_alias})) AND ({secondary_query})'
                secondary_resp_alias = _run_esearch(secondary_term_alias, remaining_secondary)
                if secondary_resp_alias:
                    pmids = secondary_resp_alias.json().get("esearchresult", {}).get("idlist", [])
                    if _collect(pmids, secondary_pmids, max_results_secondary):
                        if not alias_hits:
                            alias_hits = alias_tokens
    except Exception as e:
        print(f"ERROR: Failed to parse PubMed search results for aliases: {e}")

    # Oldest hits of the primary term. PubMed sorts dates newest first and offers
    # no ascending order, so the total is read first and the window is placed at
    # the tail. retstart is capped by the API at 9999, which only bites on
    # queries too broad to have a meaningful founding paper anyway.
    oldest_pmids = []
    try:
        if gene_clause_main and max_results_oldest > 0:
            # The phenotype clause matches exact phrases, so for a disease whose
            # name has been rewritten since the first report the primary term
            # returns nothing at all, which is precisely the case where the old
            # literature is missing. Fall back to the generic clause then.
            def _count(term, known):
                if known is not None:
                    return known
                _apply_random_delay()
                counted = _run_esearch(term, 0, sort="pub_date")
                if not counted:
                    return 0
                return int(counted.json().get("esearchresult", {}).get("count", 0))

            oldest_term = f'(({gene_clause_main})) AND ({phenotype_query})'
            total = _count(oldest_term, primary_total)
            if not total:
                oldest_term = f'(({gene_clause_main})) AND ({secondary_query})'
                total = _count(oldest_term, secondary_total)
            start = min(max(total - max_results_oldest, 0), 9999)
            if total:
                _apply_random_delay()
                tail = _run_esearch(oldest_term, max_results_oldest,
                                    sort="pub_date", retstart=start)
                if tail:
                    pmids = tail.json().get("esearchresult", {}).get("idlist", [])
                    _collect(pmids, oldest_pmids, max_results_oldest)
    except Exception as e:
        print(f"ERROR: Failed to parse PubMed oldest-first results: {e}")

    dedup_pmids = primary_pmids + secondary_pmids + oldest_pmids

    if not dedup_pmids:
        print(f"INFO: No articles found for {gene_name} and phenotypes (including aliases)")
        return [], []

    print(f"INFO: Found {len(dedup_pmids)} unique articles after alias expansion "
          f"({len(primary_pmids)} primary, {len(secondary_pmids)} secondary, "
          f"{len(oldest_pmids)} oldest)")

    summary_response = _fetch_summaries(dedup_pmids)
    if not summary_response:
        print("ERROR: Failed to fetch article summaries")
        return [], alias_hits
    
    _apply_random_delay()
    try:
        summary_data = summary_response.json()
        # Which query brought each article in. A caller can then tell an article
        # that matched the disease name from one that only matched the gene and
        # a generic word like "syndrome", which is the difference between a
        # search that worked and a search that found nothing and fell back.
        # Deliberately not `found_by`, which the prompt formatter prints.
        tiers = {}
        for pmid in oldest_pmids:
            tiers[pmid] = "oldest"
        for pmid in secondary_pmids:
            tiers[pmid] = "secondary"
        for pmid in primary_pmids:
            tiers[pmid] = "primary"
        articles = []
        for pmid in dedup_pmids:
            article_info = summary_data.get("result", {}).get(pmid, {})
            title = article_info.get("title", "No title found")
            pubdate = article_info.get("pubdate", "No date")
            articles.append({
                "pmid": pmid,
                "title": title,
                "pubdate": pubdate,
                "abstract": None,
                "query_tier": tiers.get(pmid, "")
            })
        return articles, alias_hits
    except Exception as e:
        print(f"ERROR: Failed to parse PubMed summaries: {e}")
        return [], alias_hits


def deduplicate_articles(articles):
    """
    Deduplicates a list of articles based on title. If titles are the same,
    keeps the one with the most recent publication date.
    """
    def parse_pubdate(pubdate_str):
        """Tries to parse various PubMed date formats into a datetime object."""
        formats_to_try = [
            "%Y %b %d",
            "%Y %b",
            "%Y",
        ]
        for fmt in formats_to_try:
            try:
                return datetime.strptime(pubdate_str, fmt)
            except ValueError:
                continue
        return None

    unique_articles = {}
    for article in articles:
        title = article.get("title", "").strip().lower()
        if not title:
            continue
        
        current_date = parse_pubdate(article.get("pubdate", ""))

        if title not in unique_articles:
            unique_articles[title] = article
        else:
            existing_article = unique_articles[title]
            existing_date = parse_pubdate(existing_article.get("pubdate", ""))
            
            if current_date and not existing_date:
                unique_articles[title] = article
            elif current_date and existing_date and current_date > existing_date:
                unique_articles[title] = article
    
    return list(unique_articles.values())


NOTICE_TYPES = ("Published Erratum", "Retraction of Publication",
                "Expression of Concern", "Retraction Notice")


def publication_notices(pmids):
    """The PMIDs among `pmids` that are errata or retraction notices rather
    than papers. Some journals (AJHG among them) print the erratum under the
    paper's exact title, and such a record has no abstract."""
    pmids = [str(p) for p in pmids if p]
    notices = set()
    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i + 200]
        params = {"db": "pubmed", "id": ",".join(chunk), "retmode": "json"}
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY
        resp = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", params)
        try:
            result = resp.json().get("result", {}) if resp is not None else {}
        except ValueError:
            continue
        for p in chunk:
            if any(t in NOTICE_TYPES for t in (result.get(p) or {}).get("pubtype") or ()):
                notices.add(p)
    return notices


def format_articles_for_llm(articles, use_abstracts=False, top_abstracts=None):
    """
    Format a list of articles for LLM input.
    
    Args:
        articles: List of article dictionaries
        use_abstracts: If True, includes abstract text for all articles
        top_abstracts: If set (int), only include abstracts for the first N articles.
                       This overrides use_abstracts for a hybrid mode.
                       Example: top_abstracts=10 means first 10 get abstracts, rest get titles only.
    """
    if not articles:
        return "No articles found."
    
    formatted_lines = []
    for i, article in enumerate(articles):
        title = article.get('title', 'No title')
        pmid = article.get('pmid', 'Unknown')
        found_by = article.get('found_by', [])
        abstract = article.get('abstract', '')
        
        line = f"- {title} (PMID: {pmid}) [Found by: {', '.join(found_by)}]"
        
        # Determine if we should include abstract for this article
        include_abstract = False
        if top_abstracts is not None:
            # Hybrid mode: only top N articles get abstracts
            include_abstract = (i < top_abstracts) and abstract
        elif use_abstracts:
            # All abstracts mode
            include_abstract = bool(abstract)
        
        if include_abstract:
            line += f"\n  Abstract: {abstract}"
        
        formatted_lines.append(line)
    
    return "\n".join(formatted_lines)

