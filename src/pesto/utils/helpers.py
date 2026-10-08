"""Utility functions for the application.

One function, and it matters where it sits: the package reaches load_prompt
from about forty places, so whatever this module imports is imported by every
run. It used to carry a cache reader and a numpy converter that nothing called,
and pandas came in with them, at 164 ms a process.
"""
import os

from ..config import APP_ROOT


def load_prompt(prompt_name):
    """Load a prompt from the prompts directory, by name without .txt."""
    prompt_path = os.path.join(APP_ROOT, 'prompts', f'{prompt_name}.txt')
    try:
        with open(prompt_path, 'r', encoding='utf-8') as f:
            return f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: Prompt file not found: {prompt_path}")
        return None
    except Exception as e:
        print(f"ERROR: Could not load prompt {prompt_name}: {e}")
        return None
