# PESTO

**P**rior **E**vidence **S**tatus from **T**ext and **O**ntologies. Given a gene and
a phenotype, PESTO says how much prior evidence links them, on four ordered
levels:

| Level | Meaning |
|---|---|
| Novel | no identified prior association |
| Hypothesized | mechanistic, animal or indirect support, no human genetic evidence for this pair |
| Existing | reported human genetic evidence |
| Established | replicated, accepted human genetic evidence |

Two branches answer independently, and are kept apart on purpose:

- **Literature.** PubMed is searched under the phenotype's synonyms and the gene's
  aliases, Claude Haiku 4.5 sorts every hit into four relevance bands, and Claude
  Opus 5 reads up to twenty abstracts and spreads 100 points over the four levels,
  citing the PMIDs it relied on.
- **Open Targets.** The gene's associated traits are fetched from the Open Targets
  Platform, the twenty nearest to the phenotype are ranked by an embedding model
  (BioLORD-2023), and Claude grades how each relates to the queried phenotype. The
  association score of a matching trait sets the level.

This is the code of *Large language model classifies prior evidence in
gene–phenotype associations* (Guez et al., 2026). The [`paper/`](paper/) directory
redraws every figure and table of the manuscript.

## Install

```bash
pip install pesto-genetics            # or: pipx install pesto-genetics
```

The install is about 1.3 GB, nearly all of it PyTorch, which the Open Targets
branch uses to rank traits locally. On Linux, pip's default index serves the CUDA
build of torch; nothing here uses a GPU, and the CPU build is several hundred
megabytes smaller:

```bash
pip install pesto-genetics --extra-index-url https://download.pytorch.org/whl/cpu
```

The embedding model weights (439 MB) are downloaded on first use into
`~/.cache/huggingface`. `pesto --download-models` fetches them beforehand.

## Your Anthropic API key

PESTO calls Claude through the Anthropic API, with your own key, billed to your
account. Create one at [console.anthropic.com](https://console.anthropic.com), then:

```bash
export ANTHROPIC_API_KEY=sk-ant-...
```

or put the same line, without `export`, in `~/.local/share/pesto/.env`.

**The key stays on your machine.** PESTO never writes it to disk, never prints it,
and never stores it in a saved run. It is read in one place,
[`src/pesto/credentials.py`](src/pesto/credentials.py), and handed only to the
official `anthropic` client in
[`src/pesto/services/llm_service.py`](src/pesto/services/llm_service.py), which
sends it to `api.anthropic.com` and nowhere else. Those two files are short; read
them to check.

An NCBI key (`NCBI_API_KEY`, optional) raises the PubMed rate limit from 3 to 10
requests a second and is handled the same way.

## Use

```bash
pesto --gene DDX41 --phenotype "myelodysplastic syndrome"
```

prints the literature verdict, the Open Targets verdict, how many articles were
found and read, and what the pair cost. `--json` returns everything, including
the probability distribution, the justification and the supporting PMIDs.

A table of pairs (tab-separated, columns `gene` and `phenotype`) runs in parallel:

```bash
pesto --bench pairs.tsv --workers 10      # writes pesto.tsv and pesto.cost.tsv
```

| Option | What it does |
|---|---|
| `--fulltext` | also search Europe PMC for articles that name the gene only in their body (see below) |
| `--arm NAME` | run one pipeline by name instead of choosing automatically |
| `--json` | the whole result as JSON |
| `--no-cache` | read the literature again even if this pair was already answered |
| `--resume` | with `--bench`, keep the rows an earlier attempt finished |
| `--model ID`, `--model-fast ID` | use other Claude models; the fingerprint records it |
| `--ot-encoder biolord\|sapbert\|none` | how Open Targets traits are ranked; `none` reads them all |
| `--show-pipeline` | print the blocks of the pipeline and their fingerprints, call nothing |

`pesto --help` lists the rest.

## Pipelines

A pipeline ("arm") is declared in [`src/pesto/flow/arms.py`](src/pesto/flow/arms.py)
as a list of blocks. Its fingerprint hashes every block's parameters and prompts,
and every saved run carries it, so a result says exactly which pipeline produced it.

By default PESTO asks Claude whether the phenotype is a disease, and runs:

- **`abstracts`** for a disease: PESTO as in the paper.
- **`abstracts-trait`** otherwise: the same pipeline with prompts that do not assume
  a disease, so replicated cohort and GWAS evidence can reach Established.

With `--fulltext`, it runs **`fulltext`** or **`fulltext-trait`** instead. These add a
Europe PMC search for articles that report the gene only in their body (a results
table, a supplementary file) and append the passages naming it for the reader.
PubMed cannot return those articles, and the deciding paper is often one of them.
On 20 GWAS associations whose deciding article names the gene only in its body,
`abstracts` placed 1 at Existing or above and `fulltext` 18, and both kept all 20
negative controls at Novel or Hypothesized (Supplementary Figure 3). `fulltext` was not used for the analyses of the paper,
and `fulltext-trait` has not been benchmarked.

The other arms reproduce the paper's comparisons: `titles` (titles only),
`knowledge` (the model's own knowledge, no search), and `abstracts-bare`,
`abstracts-open` and `abstracts-score` (reading-prompt ablations).

## Cost

At Anthropic list prices, the literature branch costs about $0.07–0.11 a pair
with `abstracts` (Extended Data Fig. 2 and Supplementary Figure 3), and 2–3%
more with `fulltext`; the Open Targets branch about $0.03. Well-studied genes cost more. A
pair already answered by the same pipeline is read from the cache and costs
nothing. With `--bench`, the cost of every pair is written to `pesto.cost.tsv`
as it runs.

## Where things are written

Runs are saved under `~/.local/share/pesto` (or `$PESTO_PROJECT_ROOT`), one JSON
per pair with the corpus, the prompts and the model's answer.

## Reproducing the paper

See [`paper/README.md`](paper/README.md). Every figure and table is redrawn from
the tables in `paper/results/` without calling any model. The saved model answers
behind them are on Zenodo:
[doi:10.5281/zenodo.23230222](https://doi.org/10.5281/zenodo.23230222).

## Tests

Free, no model calls, no network:

```bash
python3 tests/replay.py           # fingerprints, prompts, selection rules, parsing
python3 tests/script_imports.py   # every paper script still imports
```

## Citation

Guez J, Auwerx C, Lu W, Satterstrom FK, Berkowitz J, Fu JM, Betancur C,
Talkowski ME, Daly MJ, Karczewski KJ. Large language model classifies prior
evidence in gene–phenotype associations. 2026.
Code and data: [doi:10.5281/zenodo.23230222](https://doi.org/10.5281/zenodo.23230222).

## Licence

Code: MIT. Tables and model answers produced for the paper: CC-BY 4.0.
Third-party content they include keeps the terms of its source.
