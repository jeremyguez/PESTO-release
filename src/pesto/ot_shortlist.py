"""Ranking a gene's recorded traits by closeness to the one being asked about.

An Open Targets answer runs from a couple of dozen traits to nearly six hundred,
and a model reading all of them pays for all of them. This ranks them by cosine
similarity to the reported trait so that only the nearest few are read, which is
the whole of the saving: on the 45-pair benchmark the shortlist cuts the reading
cost by about two thirds.

Two encoders, both 12-layer and 110M parameters, both free, both CPU:

  biolord   BioLORD-2023, trained on definitions written by a language model
            over a clinical knowledge graph. Default.
  sapbert   SapBERT, trained on UMLS synonym pairs.

The choice is not cosmetic and was measured rather than assumed. On the trait
each of the 45 hand verdicts turned on, SapBERT keeps 35 of 40 inside the top 20
and BioLORD 38; both keep every literal wording of the queried trait, 25 of 25,
so the exact channel is safe either way. What SapBERT buries is the tie that
runs through a disease rather than through a word — dilated cardiomyopathy at
rank 70 for atrial fibrillation, glycogen storage disease at rank 21 for alanine
transaminase — because synonym training rewards shared wording and these share
none. Definitions mention symptoms, which is why the other model finds them.

BioLORD costs about 0.4 s more a pair, roughly 3% of the latency of the reading
call it feeds, which is why it is the default despite being the slower of the
two. Either embeds the few hundred trait names in front of it and needs nothing
but the model weights, which HuggingFace caches.
"""
from __future__ import annotations

import logging
import threading

import numpy as np

log = logging.getLogger(__name__)

BIOLORD = "biolord"
SAPBERT = "sapbert"
NONE = "none"
ENCODERS = {
    BIOLORD: "FremyCompany/BioLORD-2023",
    SAPBERT: "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",
}
# The names a caller may pass, including the one that turns the filter off.
CHOICES = (BIOLORD, SAPBERT, NONE)

_LOADED = {}
# Serialises the first load. Without it a threadpool that opens on many pairs at
# once has every worker miss the cache together and each build its own copy of a
# 439 MB model before the first finishes writing _LOADED: sixteen concurrent
# loads were seen on a twenty-worker run, a needless memory spike. With it one
# worker loads and the rest wait and reuse.
_LOAD_LOCK = threading.Lock()


def resolve_encoder(name=None):
    """The encoder to use: the argument, else the configured default."""
    from . import config
    n = (name or config.OT_ENCODER or NONE).lower()
    if n not in CHOICES:
        raise ValueError(f"unknown OT encoder {n!r}, expected one of {CHOICES}")
    return n


class Encoder:
    """One model, loaded once per process and kept.

    Pooling is a property of how the model was trained, not a preference:
    SapBERT's contrastive objective is on the [CLS] vector, while BioLORD is a
    sentence-transformers model trained on the mean of its tokens. Reading the
    wrong one costs most of the ranking quality and fails silently, which is why
    it is decided here from the model rather than left to the caller.
    """

    def __init__(self, key):
        try:
            import torch
            from transformers import AutoModel, AutoTokenizer
        except ImportError as exc:
            # Reached only when an encoder was asked for by name: without
            # torch the default is already `none`.
            raise RuntimeError(
                f"--ot-encoder {key} needs {exc.name}, which is not installed. "
                f"Install it with pip install torch transformers, or drop the "
                f"option: without them PESTO reads the whole Open Targets "
                f"answer instead (see 'Without PyTorch' in the README)."
            ) from exc
        name = ENCODERS[key]
        self.torch = torch
        self.tok = AutoTokenizer.from_pretrained(name)
        self.model = AutoModel.from_pretrained(name).eval()
        self.pool = "cls" if key == SAPBERT else "mean"

    def encode(self, texts, batch=256, max_length=32):
        out = np.zeros((len(texts), self.model.config.hidden_size),
                       dtype=np.float32)
        with self.torch.no_grad():
            for i in range(0, len(texts), batch):
                enc = self.tok(texts[i:i + batch], padding=True,
                               truncation=True, max_length=max_length,
                               return_tensors="pt")
                h = self.model(**enc).last_hidden_state
                if self.pool == "cls":
                    v = h[:, 0]
                else:
                    m = enc["attention_mask"].unsqueeze(-1).float()
                    v = (h * m).sum(1) / m.sum(1).clamp(min=1e-9)
                out[i:i + batch] = v.numpy()
        out /= np.linalg.norm(out, axis=1, keepdims=True) + 1e-12
        return out


def cached(key):
    """Whether the weights are already on this machine.

    Asks for the three files a load actually opens rather than for the whole
    snapshot: a repository holds a README and a `.gitattributes` too, and a cache
    filled by an earlier partial fetch has every weight and none of those, which
    a snapshot-level check calls missing.

    Only used to decide whether to say something before a download, so any
    failure to tell is read as "not cached" and costs one line of output.
    """
    try:
        from huggingface_hub import try_to_load_from_cache
        repo = ENCODERS[key]
        weights = any(try_to_load_from_cache(repo, f) not in (None, False)
                      for f in ("model.safetensors", "pytorch_model.bin"))
        tokens = any(try_to_load_from_cache(repo, f) not in (None, False)
                     for f in ("tokenizer.json", "vocab.txt"))
        config = try_to_load_from_cache(repo, "config.json") not in (None, False)
        return bool(weights and tokens and config)
    except Exception:
        return False


def download(names=None):
    """Fetch the weights now, rather than inside somebody's first assessment.

    The weights cannot travel in the package: 439 MB apiece against PyPI's 100 MB
    per file, and BioLORD is published under a licence of its own, derived from
    UMLS, so redistributing it is not ours to decide. Nor can a wheel run code at
    install time. So the download happens on first use, and this exists to let
    that first use be a moment somebody chose.

    Not required. Skip it and the first pair that needs a ranking pays for the
    download, once per machine, and says so while it does.
    """
    out = []
    for key in (names or (BIOLORD, SAPBERT)):
        key = key.lower()
        if key == NONE:
            continue
        if key not in ENCODERS:
            raise ValueError(f"unknown OT encoder {key!r}")
        had = cached(key)
        Encoder(key)                      # from_pretrained fetches what it needs
        out.append((key, ENCODERS[key], had))
    return out


def encoder(name=None):
    """The named encoder, built on first use and reused afterwards.

    Loading weights costs about a second, so a batch of pairs pays for it once.
    """
    key = resolve_encoder(name)
    if key == NONE:
        return None
    # Double-checked: the common path after warm-up reads _LOADED without the
    # lock; only a genuine miss takes it, and re-checks inside in case another
    # worker loaded while it waited.
    if key not in _LOADED:
        with _LOAD_LOCK:
            if key not in _LOADED:
                if not cached(key):
                    # 439 MB and twenty seconds, in the middle of what the caller
                    # asked for. Saying so beats a stalled progress bar from a
                    # library the caller never chose to use.
                    log.warning("fetching %s, 439 MB, once for this machine. "
                                "`pesto --download-models` does this beforehand.",
                                ENCODERS[key])
                _LOADED[key] = Encoder(key)
    return _LOADED[key]


def shortlist(phenotype, traits, k=None, name=None):
    """The `k` traits closest in meaning to `phenotype`, closest first.

    Returns the whole list, in the order it arrived, when the filter is off or
    when the list is already shorter than the cut. The order matters downstream:
    it is the order the reader sees, and the reader is told the entries are the
    nearest ones, so shuffling it would make that sentence false.
    """
    from . import config
    k = config.OT_TOP_K if k is None else k
    # Decided before the model is asked for, so that turning the filter off, or
    # handing it a list already shorter than the cut, does not pay a second of
    # loading weights it will never use.
    if not traits or k is None or k <= 0 or len(traits) <= k:
        return list(traits)
    enc = encoder(name)
    if enc is None:
        return list(traits)
    V = enc.encode([phenotype] + [t.get("name", "") for t in traits])
    order = np.argsort(-(V[1:] @ V[0]))[:k]
    return [traits[i] for i in order]
