"""Calls to Claude, through the Anthropic API and nothing else.

The pipeline names its models by short aliases (`haiku`, `opus5`), resolved to
API model ids in pesto.harness.MODELS. Every call lands here. The API key is
read in pesto.credentials and given to the official client, which sends it to
api.anthropic.com only.
"""
import inspect
import time
from typing import Any, Dict, Optional

import anthropic

from ..config import LLM_INITIAL_DELAY, LLM_MAX_RETRIES
from ..credentials import anthropic_api_key

VERBOSE_MODE = False

_client = None
# The 1.x SDK dropped `temperature` from messages.create, while the API still
# takes it for the models that accept one; it then goes in the request body.
_SDK_TAKES_TEMPERATURE = "temperature" in inspect.signature(
    anthropic.resources.messages.Messages.create).parameters


class LLMOverloadedError(Exception):
    """The API was still overloaded after every retry."""


class MissingAPIKeyError(RuntimeError):
    """No Anthropic API key was found in the environment or in .env."""


def _get_client():
    global _client
    if _client is None:
        key = anthropic_api_key()
        if not key:
            raise MissingAPIKeyError(
                "ANTHROPIC_API_KEY is not set. Get a key at console.anthropic.com, "
                "then: export ANTHROPIC_API_KEY=sk-ant-...")
        _client = anthropic.Anthropic(api_key=key)
    return _client


def _call_anthropic(prompt: str, model: str, temperature: float,
                    max_tokens: int = 4096) -> Dict[str, Any]:
    client_ = _get_client()
    delay = LLM_INITIAL_DELAY

    # Opus from 4.x on rejects the temperature parameter with a 400, so those
    # calls run at the provider default. Two identical Opus runs can therefore
    # answer differently; that is sampling, not drift.
    omit_temperature = model.startswith("claude-opus-")

    for attempt in range(LLM_MAX_RETRIES):
        try:
            if VERBOSE_MODE:
                print(f"INFO: Anthropic API, model {model}, "
                      f"attempt {attempt + 1}/{LLM_MAX_RETRIES}")
            kwargs = dict(model=model, max_tokens=max_tokens,
                          messages=[{"role": "user", "content": prompt}])
            if not omit_temperature:
                if _SDK_TAKES_TEMPERATURE:
                    kwargs["temperature"] = temperature
                else:
                    kwargs["extra_body"] = {"temperature": temperature}
            # The SDK refuses a long non-streamed request (a large max_tokens
            # can run past ten minutes), so those are streamed and collected.
            if max_tokens > 8192:
                with client_.messages.stream(**kwargs) as stream:
                    message = stream.get_final_message()
            else:
                message = client_.messages.create(**kwargs)
            usage = getattr(message, "usage", None)
            input_tokens = getattr(usage, "input_tokens", 0) or 0
            output_tokens = getattr(usage, "output_tokens", 0) or 0
            # Opus 5 thinks before it answers, so the first content block can be
            # a thinking block with no text. Text blocks are picked by type.
            spoken = "".join(getattr(b, "text", "") for b in message.content
                             if getattr(b, "type", "text") == "text")
            # A safety classifier can end the turn with no content, and PubMed
            # corpora reach it. Reported, so a refusal is not read as a parse bug.
            stop_reason = getattr(message, "stop_reason", None)
            if stop_reason == "refusal":
                print(f"WARNING: {model} refused to answer "
                      f"({input_tokens} input tokens)")
            return {"text": spoken, "input_tokens": input_tokens,
                    "output_tokens": output_tokens, "stop_reason": stop_reason,
                    "flex": False}
        except anthropic.APIStatusError as e:
            if e.status_code in (429, 500, 503, 529):
                if attempt < LLM_MAX_RETRIES - 1:
                    print(f"WARNING: Anthropic API returned {e.status_code}. "
                          f"Retrying in {delay}s...")
                    time.sleep(delay)
                    delay *= 2
                else:
                    raise LLMOverloadedError(
                        f"Anthropic API call failed after {LLM_MAX_RETRIES} "
                        f"attempts. Last error: {e}")
            else:
                raise

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


def call_llm_with_usage(prompt: str, model: str = "claude-haiku-4-5",
                        temperature: float = 0.0,
                        agent_name: Optional[str] = None,
                        max_tokens: Optional[int] = None) -> Dict[str, Any]:
    """One call. Returns the text and the token counts.

    `agent_name` is accepted for the callers that pass it and is not used.
    """
    result = _call_anthropic(prompt, model, temperature,
                             max_tokens=max_tokens or 4096)
    from ..cost import record
    record(model, result.get("input_tokens", 0), result.get("output_tokens", 0))
    return {"text": result.get("text"),
            "input_tokens": result.get("input_tokens", 0),
            "output_tokens": result.get("output_tokens", 0),
            "stop_reason": result.get("stop_reason"),
            "flex": False}


def call_llm(prompt: str, model: str = "claude-haiku-4-5",
             temperature: float = 0.0, agent_name: Optional[str] = None) -> str:
    """One call, returning only the text."""
    return call_llm_with_usage(prompt, model, temperature, agent_name)["text"]
