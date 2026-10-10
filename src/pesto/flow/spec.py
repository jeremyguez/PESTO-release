"""What a block is, and what a pipeline made of blocks has to satisfy.

A block declares three things about itself: the article fields it reads, the
prompts it sends, and its parameters. From those it computes a fingerprint, and
the fingerprint is the point of the exercise: a run that records the fingerprint
of every block says which pipeline produced it by reading one column.

The fingerprint hashes a block's name, parameters and prompts. It does not hash
the body of the function the block calls, in pesto.harness.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os

from ..config import APP_ROOT

PROMPTS = os.path.join(APP_ROOT, "prompts")
_SHA_CACHE = {}


def prompt_sha(name):
    """A prompt's identity. Editing the text makes old runs incomparable, which
    is a fact worth carrying rather than discovering later."""
    if name not in _SHA_CACHE:
        path = name if os.path.isabs(name) else os.path.join(PROMPTS, f"{name}.txt")
        with open(path, "rb") as fh:
            _SHA_CACHE[name] = hashlib.sha256(fh.read()).hexdigest()[:12]
    return _SHA_CACHE[name]


def digest(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:12]


def optional(default=()):
    """A parameter hashed only when set, so adding one to a block leaves every
    fingerprint computed without it unchanged."""
    return dataclasses.field(default=default, metadata={"omit_empty": True})


@dataclasses.dataclass(frozen=True)
class Block:
    """One step. Immutable, so a pipeline can be shared and never drift."""

    # Article fields this block reads. The composer checks they exist by then.
    needs = frozenset()
    # Prompt names, resolved under the package's prompts directory.
    uses = ()

    @property
    def kind(self):
        return type(self).__name__

    def params(self):
        """Everything that changes what this block does, and nothing else."""
        return {f.name: getattr(self, f.name)
                for f in dataclasses.fields(self)
                if not f.name.startswith("_")
                and not (f.metadata.get("omit_empty") and not getattr(self, f.name))}

    def takes(self, name):
        """Whether this block has a parameter of that name."""
        return any(f.name == name for f in dataclasses.fields(self))

    def fingerprint(self):
        return digest({"kind": self.kind, "params": self.params(),
                       "prompts": {p: prompt_sha(p) for p in self.uses}})

    def but(self, **changes):
        return dataclasses.replace(self, **changes)

    def describe(self):
        shown = ", ".join(f"{k}={v!r}" for k, v in sorted(self.params().items())
                          if v is not None and v != ())
        return f"{self.kind}({shown})"


@dataclasses.dataclass(frozen=True)
class Reader(Block):
    """Turns a corpus into evidence. Declares which kind, so a judge can be
    matched to it before anything is paid for."""
    evidence = "none"


@dataclasses.dataclass(frozen=True)
class Judge(Block):
    """Turns evidence into a verdict, by rule alone. No model, no network, so a
    judgement can be replayed from disk and argued with."""
    reads = "none"

    def garbled(self, text):
        """Whether the reader's answer is too malformed to be trusted, in which
        case the run asks the reader again."""
        return False


@dataclasses.dataclass(frozen=True)
class Arm:
    """A pipeline, written down rather than spelled out in a script.

    Architecture is what is in these lists; parameters are the numbers inside
    the blocks. Changing `sources` is a different pipeline. Changing a quota
    inside one of them is the same pipeline, differently set, and the
    fingerprint tells the two apart without anyone having to remember which.
    """
    name: str
    expand: tuple
    sources: tuple
    dedupe: Block
    hydrate: Block
    annotate: tuple
    select: Block
    read: Reader
    judge: Judge
    # Where the abstracts are fetched, counted in annotators already run. Zero
    # is before all of them; one places it after the title pass, which fetches
    # only for what survived.
    hydrate_at: int = 0
    note: str = ""
    # The last publication date searched, YYYY/MM/DD; empty searches everything.
    max_date: str = ""

    def blocks(self):
        annotate = list(self.annotate)
        annotate.insert(min(self.hydrate_at, len(annotate)), self.hydrate)
        return (tuple(self.expand) + tuple(self.sources) + (self.dedupe,)
                + tuple(annotate) + (self.select, self.read, self.judge))

    def fingerprint(self):
        # Hashed only when set, so every run made without a date keeps its
        # fingerprint and stays in the cache.
        parts = [b.fingerprint() for b in self.blocks()]
        if self.max_date:
            parts.append({"max_date": self.max_date})
        return digest(parts)

    def but(self, name=None, **changes):
        return dataclasses.replace(self, name=name or self.name, **changes)

    def using(self, reader=None, worker=None, max_date=None, notes=()):
        """The same arm read by a different model, over an older literature,
        or with notes added to the reading prompt.

        Swapping models is a change of parameter and not of architecture, so
        the arm stays the arm, and the fingerprint moves, which is the honest
        record that the numbers may too. A date ceiling and notes move it the
        same way.
        """
        def swap(block, model):
            # A block whose model is None is not calling anything: that is how
            # a block says it replays a recorded synonym list. Giving it a model
            # would make it expand live, which is a different arm.
            return (block.but(model=model)
                    if model and block.params().get("model") else block)

        read = swap(self.read, reader or worker)
        if notes:
            if not read.takes("notes"):
                raise ValueError(f"{self.name}: {read.kind} reads no articles, "
                                 f"so its prompt has no place for a note")
            read = read.but(notes=tuple(read.notes) + tuple(notes))
        return dataclasses.replace(
            self,
            expand=tuple(swap(b, worker) for b in self.expand),
            annotate=tuple(swap(b, worker) for b in self.annotate),
            read=read,
            max_date=max_date or self.max_date)

    def validate(self):
        """Everything checkable before a single call is made.

        Two failures are worth catching here rather than in the middle of a run
        that costs sixteen dollars: a judge that cannot read what its reader
        produces, and a block that wants abstracts before they are fetched.
        """
        if self.judge.reads != self.read.evidence:
            raise ValueError(
                f"{self.name}: {self.judge.kind} reads {self.judge.reads!r} but "
                f"{self.read.kind} produces {self.read.evidence!r}")
        from .types import SEARCHED
        have = set(SEARCHED)
        for i, block in enumerate(tuple(self.annotate) + (self.select, self.read)):
            if i >= self.hydrate_at:
                have |= {"abstract"}
            missing = set(block.needs) - have
            if missing:
                raise ValueError(
                    f"{self.name}: {block.kind} reads {sorted(missing)}, which "
                    f"nothing upstream provides at that point")
        return self

    def describe(self):
        lines = [f"{self.name}  [{self.fingerprint()}]"]
        if self.max_date:
            lines.append(f"  literature published up to {self.max_date}")
        for block in self.blocks():
            lines.append(f"  {block.fingerprint()}  {block.describe()}")
        return "\n".join(lines)
