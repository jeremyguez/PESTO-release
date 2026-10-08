"""Gene names: the HGNC aliases and approved name a gene is published under.

Read from aliases.tsv, shipped with the package, and built once on first use.
"""
import os
import re

import pandas as pd

from ..config import APP_ROOT

ALIASES_FILE = os.path.join(APP_ROOT, "aliases.tsv")


class DataLoader:
    """The alias table, loaded lazily."""

    def __init__(self):
        self.alias_map = None  # Maps uppercase alias/previous symbol -> approved symbol
        self.approved_to_aliases = None  # Maps approved symbol uppercase -> set of tokens (including approved)
        self.approved_names = None  # Maps approved symbol uppercase -> HGNC approved name


    def _load_alias_map(self):
        """Lazy-load HGNC aliases TSV and build alias->approved map (first hit wins)."""
        if self.alias_map is not None:
            return

        if not os.path.exists(ALIASES_FILE):
            print(f"WARNING: aliases file not found at {ALIASES_FILE}")
            self.alias_map = {}
            return

        try:
            df = pd.read_csv(ALIASES_FILE, sep='\t', dtype=str).fillna('')
        except Exception as e:
            print(f"ERROR: Could not load aliases file {ALIASES_FILE}: {e}")
            self.alias_map = {}
            return

        alias_map = {}
        approved_to_aliases = {}
        approved_names = {}

        def add_alias(token, approved_symbol):
            alias_upper = token.strip().upper()
            if not alias_upper:
                return
            if alias_upper not in alias_map:
                alias_map[alias_upper] = approved_symbol
            approved_upper = approved_symbol.upper()
            approved_to_aliases.setdefault(approved_upper, set()).add(alias_upper)

        for _, row in df.iterrows():
            approved = (row.get('Approved symbol') or '').strip()
            if not approved:
                continue
            approved_symbol = approved
            add_alias(approved_symbol, approved_symbol)

            name = (row.get('Approved name') or '').strip()
            if name:
                approved_names[approved_symbol.upper()] = name

            prev_symbols = (row.get('Previous symbols') or '')
            alias_symbols = (row.get('Alias symbols') or '')

            for cell in (prev_symbols, alias_symbols):
                if not cell:
                    continue
                for token in re.split(r'[,;]', cell):
                    add_alias(token, approved_symbol)

        self.alias_map = alias_map
        self.approved_to_aliases = approved_to_aliases
        self.approved_names = approved_names

    def _resolve_gene_symbol(self, gene_name):
        """Return approved symbol if gene_name matches an alias/previous symbol."""
        if not gene_name:
            return gene_name
        self._load_alias_map()
        return self.alias_map.get(gene_name.strip().upper(), gene_name)

    def get_aliases_for_gene(self, gene_name):
        """
        Return list of aliases/previous symbols for the approved gene (excludes approved symbol).
        """
        if not gene_name:
            return []
        self._load_alias_map()
        resolved = self._resolve_gene_symbol(gene_name)
        approved_upper = (resolved or "").strip().upper()
        if not approved_upper:
            return []
        tokens = list(self.approved_to_aliases.get(approved_upper, set()))
        # Remove the approved symbol itself (keep only alternatives)
        aliases = [t for t in tokens if t != approved_upper]
        return aliases

    def get_search_names_for_gene(self, gene_name):
        """Alias symbols plus the HGNC approved name, for literature search.

        Authors of the papers a curator cites often never write the symbol: the
        dystrophin literature says dystrophin, and of the eighteen PMIDs GenCC
        cites for DMD, five carry the symbol and fifteen carry the name. The name
        is a phrase rather than a token, so it belongs with the aliases only when
        searching text, not when resolving a symbol to its approved form.
        """
        if not gene_name:
            return []
        aliases = self.get_aliases_for_gene(gene_name)
        approved = (self._resolve_gene_symbol(gene_name) or "").strip().upper()
        name = (self.approved_names or {}).get(approved, "")
        # A name that merely spells the symbol out adds a synonym of itself and
        # doubles the query for nothing.
        if name and name.strip().upper() not in {approved, *(a.upper() for a in aliases)}:
            return aliases + [name]
        return aliases


data_loader = DataLoader()
