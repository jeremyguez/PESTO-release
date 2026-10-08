"""The API keys: where they are read, and what is done with them.

PESTO needs one secret, an Anthropic API key, and can use a second, an NCBI key
that raises the PubMed rate limit from 3 to 10 requests a second. Both are read
here and nowhere else, from either of:

  - the environment variables ANTHROPIC_API_KEY and NCBI_API_KEY;
  - a .env file you wrote yourself in the PESTO data directory
    (~/.local/share/pesto/.env, or $PESTO_PROJECT_ROOT/.env when that is set).
    A variable already set in the shell wins over the file.

PESTO never writes a key to disk, never prints it, and never records it in a
saved run. The Anthropic key is handed to the official `anthropic` client in
services/llm_service.py, which sends it to api.anthropic.com and nowhere else.
The NCBI key is sent to eutils.ncbi.nlm.nih.gov as the `api_key` parameter NCBI
asks for, in services/pubmed_service.py. Those three files are the whole audit.
"""
import os

from dotenv import load_dotenv

_loaded = False


def _load():
    global _loaded
    if not _loaded:
        from .config import PROJECT_ROOT
        load_dotenv(os.path.join(PROJECT_ROOT, ".env"), override=False)
        _loaded = True


def anthropic_api_key():
    _load()
    return os.environ.get("ANTHROPIC_API_KEY") or None


def ncbi_api_key():
    _load()
    return os.environ.get("NCBI_API_KEY") or None
