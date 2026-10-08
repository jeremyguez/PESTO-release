"""Configuration for the novelty-triage pipeline.

Two roots, kept apart on purpose. APP_ROOT is the code as installed: prompts,
the alias table, everything shipped alongside the modules. PROJECT_ROOT is where
a run reads its inputs and writes its outputs, which stops being the same place
the moment the package is installed into site-packages rather than run from a
checkout. It is read from PESTO_PROJECT_ROOT, and otherwise falls back to a
per-user data directory, so a fresh install writes somewhere sane instead of
into the interpreter. The scripts behind the paper set that variable to the
paper/ directory, which is what keeps their paths unchanged.

API keys are not read here: see credentials.py.
"""
import importlib.util
import os

# --- Path Configuration ---
APP_ROOT = os.path.dirname(os.path.realpath(__file__))


def _default_project_root():
    xdg = os.environ.get("XDG_DATA_HOME") or \
        os.path.join(os.path.expanduser("~"), ".local", "share")
    return os.path.join(xdg, "pesto")


PROJECT_ROOT = os.path.abspath(
    os.environ.get("PESTO_PROJECT_ROOT") or _default_project_root())

# --- URLs ---
BASE_URL_NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# --- File Paths ---
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RAW_DIR = os.path.join(DATA_DIR, 'raw')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')

# Where a single pair's evidence is written.
RUNS_BASE_DIR = os.path.join(DATA_DIR, 'runs')
NOVEL_RUNS_DIR = os.path.join(RUNS_BASE_DIR, 'novel')

# Reference tables read by the paper scripts.
GENCC_FILE = os.path.join(RAW_DIR, 'gencc-submissions.tsv')
AOU_FILE = os.path.join(RAW_DIR, 'AoU_results.tsv')
BRAVA_FILE = os.path.join(RAW_DIR, 'Duncan_results.tsv')
ASD_REVIEW_XLSX = os.path.join(RAW_DIR, 'ASD_manual_review_CB_CA_20260406.xlsx')
ASD_RESULTS_TSV = os.path.join(RAW_DIR, 'ASD_novel_associations_results.tsv')

# --- Run parameters ---
# The cheap model, used where no arm says otherwise (the Open Targets sieve).
DEFAULT_MODEL = "claude-haiku-4-5"
DEFAULT_TEMPERATURE = 0.0
DEFAULT_MAX_PRIMARY = 30
DEFAULT_MAX_SECONDARY = 10
# Zero: the NCBI rate limit is enforced once, in pubmed_service, for the
# whole client.
DEFAULT_DELAY_MS = 0
DEFAULT_WORKERS = 20
DEFAULT_SEED = 42

# --- How the Open Targets answer is searched for the reported trait ---
#
# "tagged_shortlist" ranks the answer by an embedding model, sieves it with
# claude-haiku and grades the survivors one by one with claude-opus. See
# ot_matcher.py.
OT_MATCHER = "tagged_shortlist"
# Which embedding model ranks the answer: "biolord", "sapbert", or "none" to
# read the whole list. See ot_shortlist.py for what the choice costs and buys.
# An install without PyTorch cannot rank, so it reads the whole list instead;
# the encoder is part of the cache key, so those runs are never mistaken for
# ranked ones.
TORCH_INSTALLED = (importlib.util.find_spec("torch") is not None and
                   importlib.util.find_spec("transformers") is not None)
OT_ENCODER = os.environ.get("OT_ENCODER", "biolord" if TORCH_INSTALLED else "none")
# How many of the nearest traits are read. Zero or less reads the whole answer.
OT_TOP_K = int(os.environ.get("OT_TOP_K", "20"))
# Whether claude-haiku sieves the shortlist before claude-opus grades it. It
# removes about four traits of twenty for a twentieth of a cent, and dropped
# nothing that carried a tie across 45 benchmark pairs.
OT_GATE = os.environ.get("OT_GATE", "1") not in ("0", "false", "False", "")
# Whether each tag is accompanied by a reason. Off by default: it is 60% of the
# cost of the reading call. Turn it on when the grades themselves are read.
OT_REASONS = os.environ.get("OT_REASONS", "0") not in ("0", "false", "False", "")
# The tagged-match prompt that grades the shortlist. The short name enters the
# cache key, so a run under one wording is never reused for another.
OT_MATCH_PROMPTS = {"v7": "ot_tagged_match_v7_prompt"}
OT_MATCH_PROMPT = "v7"


def ot_match_prompt_name(reasons=None, prompt=None):
    """The prompt file the tagged reader is sent."""
    if OT_REASONS if reasons is None else reasons:
        return "ot_tagged_match_reasons_prompt"
    which = OT_MATCH_PROMPT if prompt is None else prompt
    return OT_MATCH_PROMPTS.get(which, OT_MATCH_PROMPTS["v7"])


def ot_agent_version(matcher=None, encoder=None, top_k=None, reasons=None,
                     prompt=None):
    """What the Open Targets agents saw, excluding verdict rules.

    TIES, score bands and whether the gate is applied afterwards can be
    replayed from a saved record. The encoder, the cut, the reasons flag
    and which match prompt was sent change which traits were graded, so
    those belong here.
    """
    m = matcher or OT_MATCHER
    k = OT_TOP_K if top_k is None else top_k
    which = OT_MATCH_PROMPT if prompt is None else prompt
    return "%s:%s:%d%s:%s" % (
        m, (encoder or OT_ENCODER), max(k, 0),
        ":reasons" if (OT_REASONS if reasons is None else reasons) else "",
        which)


def ot_match_version(matcher=None, encoder=None, top_k=None, gate=None,
                     reasons=None):
    """What makes two Open Targets answers comparable.

    A cached answer may be reused only if it was produced under the same
    settings: the encoder, the cut and the sieve all change which traits were
    read and therefore what the answer says.
    """
    return ot_agent_version(matcher, encoder, top_k, reasons) + (
        ":gate" if (OT_GATE if gate is None else gate) else "")


OT_MATCH_VERSION = ot_match_version()

# --- Cost accounting ---
# Anthropic list prices, USD per million tokens (input, output).
PRICE_PER_MTOK_INPUT = 1.00
PRICE_PER_MTOK_OUTPUT = 5.00
MODEL_PRICES = {
    "haiku": (1.00, 5.00),
    "claude-haiku-4-5": (1.00, 5.00),
    "opus": (5.00, 25.00),
    "opus5": (5.00, 25.00),
    "claude-opus-5": (5.00, 25.00),
}


def price_usd(model, input_tokens, output_tokens):
    """List price of one call. Unknown models are charged as Haiku."""
    key = (model or "").split("/")[-1].split("@")[0]
    pin, pout = MODEL_PRICES.get(key, (PRICE_PER_MTOK_INPUT, PRICE_PER_MTOK_OUTPUT))
    return (int(input_tokens or 0) / 1e6) * pin + (int(output_tokens or 0) / 1e6) * pout

# --- API Settings ---
LLM_MAX_RETRIES = 3
LLM_INITIAL_DELAY = 2
API_REQUEST_TIMEOUT = 30

VERDICT_ORDER = ["Novel", "Hypothesized", "Existing", "Established"]
VERDICT_COLORS = {
    "Novel": "#b91c1c",
    "Hypothesized": "#b45309",
    "Existing": "#1d4ed8",
    "Established": "#15803d",
}
