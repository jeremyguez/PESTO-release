"""Parsers for LLM responses: the synonym list, and a novelty assessment."""
import re


def parse_phenotype_synonyms_response(response_text):
    """Parses the XML response from the phenotype synonyms LLM."""
    if not response_text:
        return []
    
    try:
        # Try to extract synonyms using regex
        synonyms_match = re.search(r'<synonyms>(.*?)</synonyms>', response_text, re.DOTALL | re.IGNORECASE)
        
        if not synonyms_match:
            print("WARNING: Could not find <synonyms> tags in response")
            return []
        
        synonyms_content = synonyms_match.group(1).strip()
        
        if not synonyms_content:
            print("INFO: No synonyms returned (empty content)")
            return []
        
        # Extract individual synonym tags
        synonym_matches = re.findall(r'<synonym>(.*?)</synonym>', synonyms_content, re.DOTALL | re.IGNORECASE)
        
        synonyms = [s.strip() for s in synonym_matches if s.strip()]
        
        print(f"INFO: Parsed {len(synonyms)} phenotype synonyms")
        return synonyms
        
    except Exception as e:
        print(f"ERROR: Failed to parse phenotype synonyms response: {e}")
        return []


def parse_novelty_assessment_response(response_text, all_available_pmids):
    """Parses the XML response from the novelty assessment LLM."""
    if not response_text:
        return None
    
    try:
        # Extract verdict
        verdict_match = re.search(r'<verdict>(.*?)</verdict>', response_text, re.DOTALL | re.IGNORECASE)
        verdict = verdict_match.group(1).strip() if verdict_match else "Unknown"
        
        # Validate verdict
        valid_verdicts = ["Novel", "Hypothesized", "Existing", "Established"]
        if verdict not in valid_verdicts:
            print(f"WARNING: Invalid verdict '{verdict}', defaulting to 'Unknown'")
            verdict = "Unknown"
        
        # Extract justification
        just_match = re.search(r'<justification>(.*?)</justification>', response_text, re.DOTALL | re.IGNORECASE)
        justification = just_match.group(1).strip() if just_match else "No justification provided."
        
        # Extract PMIDs
        pmids_match = re.search(r'<pmids>(.*?)</pmids>', response_text, re.DOTALL | re.IGNORECASE)
        pmids = []
        
        if pmids_match:
            pmids_content = pmids_match.group(1).strip()
            # Extract individual PMID tags
            pmid_matches = re.findall(r'<pmid>(.*?)</pmid>', pmids_content, re.DOTALL | re.IGNORECASE)
            
            for pmid in pmid_matches:
                pmid = pmid.strip()
                # Extract digits only
                digit_match = re.search(r'\d+', pmid)
                if digit_match:
                    extracted_pmid = digit_match.group(0)
                    # Validate against available PMIDs
                    if extracted_pmid in all_available_pmids:
                        pmids.append(extracted_pmid)
                    else:
                        print(f"WARNING: PMID {extracted_pmid} not in available PMIDs list")
        
        # Limit to 10 PMIDs
        pmids = pmids[:10]
        
        print(f"INFO: Parsed novelty assessment - Verdict: {verdict}, PMIDs: {len(pmids)}")
        
        return {
            "verdict": verdict,
            "justification": justification,
            "pmids": pmids
        }
        
    except Exception as e:
        print(f"ERROR: Failed to parse novelty assessment response: {e}")
        return None


# ---------------------------------------------------------------------------
# Probabilistic mode parsers
# ---------------------------------------------------------------------------

PENETRANCE_PROBA_KEYS = ["mendelian", "high", "moderate", "complex"]
INHERITANCE_PROBA_KEYS = ["dominant", "inc_dom", "incomplete", "codominant", "inc_rec", "recessive"]
ONSET_PROBA_KEYS = ["prenatal", "neonatal", "infancy", "childhood", "adolescence", "adulthood", "late"]
SEVERITY_PROBA_KEYS = ["lethal", "severe", "moderate", "mild", "verymild"]


# ---------------------------------------------------------------------------
# Phenotyping parsers
# ---------------------------------------------------------------------------


