"""Selected GO + PEPPER bars for figure 3c (ASD) and 3d (SCZ).

Universe is the 450 genes scored on both phenotypes. FDR is BH over
all GO terms (or all continuous PEPPER columns, embeddings included).
The figure keeps a non-redundant subset; everything else at FDR < 0.1
is written to results/fig3_cd_leftovers.tsv.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
# PEPPER gene features (Guez et al. 2026) are not redistributed here; see
# scripts/92_fig3_random50.py.
PEPPER = ROOT / "data" / "pepper"
GMT = {
    "BP": ROOT / "data" / "go" / "GO_Biological_Process_2026.gmt",
    "MF": ROOT / "data" / "go" / "GO_Molecular_Function_2026.gmt",
    "CC": ROOT / "data" / "go" / "GO_Cellular_Component_2026.gmt",
    "SynGO": ROOT / "data" / "syngo" / "SynGO_2024.gmt",
}
ID = {"ensg", "hgnc", "chrom", "gene_symbol", "gene"}
PCA = ("human_", "mouse_", "protdim_", "coexp_")

# Least-redundant, most-significant representatives. Order is display order
# (most significant first within kind).
ASD = [
    ("GO", "Regulation of DNA-templated Transcription (GO:0006355)",
     "DNA-templated transcription"),
    ("GO", "Chromatin Remodeling (GO:0006338)",
     "chromatin remodeling"),
    ("GO", "Negative Regulation of Gene Expression, Epigenetic (GO:0045814)",
     "epigenetic gene silencing"),
    ("GO", "Nervous System Development (GO:0007399)",
     "nervous system development"),
    ("GO", "Integral Component Of Postsynaptic Membrane (GO:0099055) CC",
     "postsynaptic membrane"),
    ("PEPPER", "phastCons17way_max", "conservation (phastCons)"),
    ("PEPPER", "transcript_length", "transcript length"),
    ("PEPPER", "5UTR_length", "5' UTR length"),
    ("PEPPER", "promoter_count", "FANTOM promoter count"),
]
SCZ = [
    ("GO", "Postsynaptic Density (GO:0014069)",
     "postsynaptic density"),
    ("GO", "Integral Component Of Presynaptic Membrane (GO:0099056) CC",
     "presynaptic membrane"),
    ("GO", "Glutamate Receptor Signaling Pathway (GO:0007215)",
     "glutamate receptor signaling"),
    ("GO", "Modulation of Chemical Synaptic Transmission (GO:0050804)",
     "modulation of transmission"),
    ("GO", "Regulation of Neuronal Synaptic Plasticity (GO:0048168)",
     "neuronal synaptic plasticity"),
    ("GO", "Learning (GO:0007612)", "learning"),
    ("PEPPER", "phastCons17way_max", "conservation (phastCons)"),
    ("PEPPER", "transcript_length", "transcript length"),
    ("PEPPER", "5UTR_length", "5' UTR length"),
    ("PEPPER", "promoter_count", "FANTOM promoter count"),
]
# Bipolar: glutamate plus the excitability and morphology terms that
# distinguish it from the SCZ synapse set, then the PEPPER features that
# actually move (transcript architecture, tissue specificity, PPI degree).
BIP = [
    ("GO", "Glutamate Receptor Signaling Pathway (GO:0007215)",
     "glutamate receptor signaling"),
    ("GO", "Dendrite (GO:0030425)",
     "dendrite"),
    ("GO", "Monoatomic Cation Transmembrane Transport (GO:0098655)",
     "cation transmembrane transport"),
    ("GO", "Action Potential (GO:0001508)",
     "action potential"),
    ("GO", "Synapse Assembly (GO:0007416)",
     "synapse assembly"),
    ("GO", "Axon Guidance (GO:0007411)",
     "axon guidance"),
    ("PEPPER", "transcript_length", "transcript length"),
    ("PEPPER", "5UTR_length", "5' UTR length"),
    ("PEPPER", "tau", "tissue specificity (tau)"),
    ("PEPPER", "PPI_degree_quantile", "PPI degree"),
]


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=0):
        k = n - rank
        prev = min(prev, pvals[i] * n / k)
        q[i] = min(1.0, prev)
    return q


def load_pheno(pheno):
    by = {}
    for path in (
        ROOT / "benchmark" / "asc_ref299" / "pesto.tsv",
        ROOT / "benchmark" / "asd_extra101" / "pesto.tsv",
        ROOT / "benchmark" / "scz_extra101" / "pesto.tsv",
        ROOT / "benchmark" / "extra50_random" / "pesto.tsv",
    ):
        with path.open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("phenotype") == pheno:
                    by[r["gene"].upper()] = r["lit_mean"]
    return by


def reaching(by):
    universe = set(by)
    fg = {g for g in universe if by[g] in ("Existing", "Established")}
    return fg, universe - fg, universe


def read_gmt(universe):
    terms = []
    for path in GMT.values():
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                genes = {g.upper() for g in f[2:] if g} & universe
                if len(genes) >= 5:
                    terms.append((f[0], genes))
    return terms


def go_rows(fg, universe):
    rows = []
    for name, genes in read_gmt(universe):
        a = len(fg & genes)
        b = len(fg) - a
        c = len((universe - fg) & genes)
        d = len(universe - fg) - c
        if a == 0:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        fold = (a / len(fg)) / ((a + c) / len(universe))
        rows.append({
            "kind": "GO", "item": name, "k": a, "K": a + c,
            "n_fg": len(fg), "N": len(universe), "effect": fold, "p": p,
        })
    qs = bh([r["p"] for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = q
    return rows


def classify(s):
    v = s.dropna()
    if v.empty:
        return None
    if set(np.unique(v.to_numpy())) <= {0, 1, 0.0, 1.0}:
        return "binary"
    if pd.api.types.is_numeric_dtype(s):
        return "continuous"
    return None


def pepper_rows(fg, universe, feat, mapping):
    m = mapping[mapping.gene.isin(universe)].drop_duplicates("ensg")
    tab = feat.merge(m[["ensg", "gene"]], on="ensg").drop_duplicates("gene")
    rows = []
    for col in tab.columns:
        if col in ID:
            continue
        if classify(tab[col]) != "continuous":
            continue
        vals = {g: v for g, v in zip(tab.gene, tab[col]) if pd.notna(v)}
        a = [vals[g] for g in fg if g in vals]
        b = [vals[g] for g in (universe - fg) if g in vals]
        if len(a) < 8 or len(b) < 8:
            continue
        U, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = (2 * U) / (len(a) * len(b)) - 1
        rows.append({
            "kind": "PEPPER", "item": col, "k": np.nan, "K": np.nan,
            "n_fg": len(a), "N": len(a) + len(b), "effect": rb, "p": p,
        })
    qs = bh([r["p"] for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = q
    return rows


def pack(contrast, panel, selected, go, pep):
    by = {r["item"]: r for r in go + pep}
    bars, leftovers = [], []
    keep = {(kind, item) for kind, item, _ in selected}
    for i, (kind, item, label) in enumerate(selected):
        r = by[item]
        bars.append({
            "panel": panel, "contrast": contrast, "kind": kind,
            "item": item, "label": f"{kind} · {label}",
            "order": i, "q": r["q"], "p": r["p"],
            "effect": r["effect"], "nlog10q": -np.log10(max(r["q"], 1e-300)),
            "k": r["k"], "K": r["K"], "n_fg": r["n_fg"], "N": r["N"],
        })
    for r in go + pep:
        if r["q"] >= 0.1:
            continue
        if (r["kind"], r["item"]) in keep:
            continue
        if r["kind"] == "PEPPER" and r["item"].startswith(PCA):
            continue
        leftovers.append({
            "panel": panel, "contrast": contrast, "kind": r["kind"],
            "item": r["item"], "q": r["q"], "p": r["p"],
            "effect": r["effect"], "k": r["k"], "K": r["K"],
            "n_fg": r["n_fg"], "N": r["N"],
        })
    leftovers.sort(key=lambda r: (r["kind"], r["q"]))
    return bars, leftovers


def main():
    mapping = pd.read_csv(
        PEPPER / "obs_exp_for_loeuf_missense.tsv",
        sep="\t", usecols=["gene_symbol", "ensg"],
    )
    mapping = mapping[mapping.gene_symbol.notna() & mapping.ensg.notna()]
    mapping["gene"] = mapping.gene_symbol.str.upper()
    feat = pd.read_csv(
        PEPPER / "gene_features_for_s_het.tsv.gz",
        sep="\t", compression="gzip",
    )

    bars, leftovers = [], []
    for pheno, short, selected, panel in (
        ("Autism Spectrum Disorder", "ASD", ASD, "c"),
        ("schizophrenia", "SCZ", SCZ, "d"),
    ):
        fg, rest, universe = reaching(load_pheno(pheno))
        print(f"{short}: reaching {len(fg)}  rest {len(rest)}  N={len(universe)}")
        go = go_rows(fg, universe)
        pep = pepper_rows(fg, universe, feat, mapping)
        b, left = pack(
            f"{short}-reaching vs rest ({len(fg)} vs {len(rest)})",
            panel, selected, go, pep,
        )
        bars += b
        leftovers += left
        print(f"  GO FDR<0.1 {sum(r['q']<0.1 for r in go)}  "
              f"PEPPER FDR<0.1 {sum(r['q']<0.1 for r in pep)}  "
              f"shown {len(b)}  leftover non-embedding {len(left)}")

    dest = ROOT / "results" / "fig3_cd_bars.tsv"
    pd.DataFrame(bars).to_csv(dest, sep="\t", index=False)
    left_dest = ROOT / "results" / "fig3_cd_leftovers.tsv"
    pd.DataFrame(leftovers).to_csv(left_dest, sep="\t", index=False)
    print(f"wrote {dest}")
    print(f"wrote {left_dest}")


if __name__ == "__main__":
    main()
