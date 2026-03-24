# reconparser

**reconparser** is a Python library for parsing phylogenetic reconciliation outputs from various software tools.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **ALE (Amalgamated Likelihood Estimation)** — Parser for classic ALE reconciliation outputs
  - Consensus gene trees (`.ucons_tree`)
  - Transfer events with frequencies (`.uTs`)
  - Maximum likelihood reconciliations (`.uml_rec`)
  - ML rates (Duplications, Transfers, Losses)
  - Log-likelihood values
  - Per-branch reconciliation statistics
  - Support for ALE v0.4 and v1.0 output formats

- **AleRax** — Parser for AleRax reconciliation outputs (v1.2+)
  - Run-level access (`AleRaxRun`): species trees, model parameters, likelihoods, global transfers, per-species events, origination probabilities, coverage data
  - Per-family access (`AleRaxFamily`): consensus trees, sampled gene trees, annotated reconciled trees, transfers, event counts, per-species events
  - Fully lazy: nothing is loaded until you ask for it, making it efficient for large datasets with thousands of species

- **Extensible** — Designed to support additional reconciliation tools in future releases

## Installation

```bash
pip install reconparser  # (not available on PyPI yet)
```

Or install from source:

```bash
git clone https://github.com/aadavin/reconparser.git
cd reconparser
pip install -e .
```

### Dependencies

- Python ≥ 3.10
- ete3 (for tree handling)
- pandas (for data tables)
- numpy

## Quick Start

### Parsing ALE output

```python
from reconparser import ALEParser

parser = ALEParser("results.ale")

# Consensus gene tree
gene_tree = parser.get_consensus_tree()
print(f"Gene tree has {len(gene_tree.get_leaves())} leaves")

# Transfer events
transfers = parser.get_transfers()
print(f"Found {len(transfers)} transfer events")
print(transfers.head())

# ML reconciliation rates
ml_rates = parser.get_ml_rates()
print(f"D: {ml_rates['duplications']:.4f}  "
      f"T: {ml_rates['transfers']:.4f}  "
      f"L: {ml_rates['losses']:.4f}")

# Log-likelihood
print(f"Log-likelihood: {parser.get_log_likelihood():.2f}")

# Per-branch statistics
branch_stats = parser.get_branch_statistics()
print(branch_stats[branch_stats['transfers'] > 0])

# Reconciled gene trees (default 100)
gene_trees = parser.get_reconciled_gene_trees()
print(f"Parsed {len(gene_trees)} reconciled gene trees")
```

### Parsing AleRax output

AleRax output is accessed at two levels: the full run and individual gene families.

```python
from reconparser import AleRaxRun

# Point to the AleRax output directory
run = AleRaxRun("path/to/alerax_output")

# Run metadata
info = run.get_run_info()
print(f"AleRax {info['version']}  |  "
      f"{info['num_families']} families  |  "
      f"{info['num_species']} species")

# Species tree
species_tree = run.get_species_tree()

# Model parameters (DTL rates per gene)
params = run.get_model_parameters()

# Per-family likelihoods
lk = run.get_per_family_likelihoods()
print(f"Total log-likelihood: {run.get_total_log_likelihood():.2f}")

# Global transfers across all families
transfers = run.get_transfers()
print(transfers.head(5))
```

Drilling into a single gene family:

```python
fam = run.get_family("K00192")

# Consensus gene tree (with support values)
consensus = fam.get_consensus_tree()

# Averaged transfers across all samples
fam_transfers = fam.get_transfers()

# Event counts for a single sample (lazy — only reads that file)
counts = fam.get_event_counts(sample=0)
# {'S': 63, 'SL': 65, 'D': 5, 'DL': 0, 'T': 33, 'TL': 14, 'L': 0, 'Leaf': 102}

# Event counts across all samples as a DataFrame
all_counts = fam.get_all_event_counts()
print(all_counts[["S", "D", "T"]].describe())

# Load a single sampled gene tree without loading all of them
tree_5 = fam.get_sampled_gene_tree(5)

# Annotated reconciled gene tree (.rec_uml format)
rec_tree = fam.get_reconciled_gene_tree()
# Node names contain .T@donor->recipient, .D@species, .S annotations
```

## Detailed Usage

### ALE output files

ALE produces three main output files:

1. **`.ucons_tree`** — Consensus gene tree in Newick format
2. **`.uTs`** — Transfer events (tab-separated: from, to, frequency)
3. **`.uml_rec`** — Maximum likelihood reconciliation containing the species tree, reconciled gene trees (variable number, default 100), ML rate estimates (D, T, L), log-likelihood, summary statistics, and per-branch event counts

#### Analyzing ALE transfers

```python
from reconparser import ALEParser

parser = ALEParser("my_analysis.ale")
transfers = parser.get_transfers()

# Filter high-frequency transfers
high_freq = transfers[transfers['freq'] > 0.1]

# Top 10 by frequency
for _, row in transfers.nlargest(10, 'freq').iterrows():
    print(f"{row['from']} -> {row['to']}: {row['freq']:.3f}")
```

#### ALE branch-level analysis

```python
parser = ALEParser("my_analysis.ale")
branch_stats = parser.get_branch_statistics()

# Find hotspot branches
hotspots = branch_stats[
    (branch_stats['transfers'] > 2) | (branch_stats['duplications'] > 2)
]
print(hotspots[['branch_id', 'duplications', 'transfers', 'losses']])
```

### AleRax output directory

AleRax produces a directory structure containing run-level and per-family files:

```
output_dir/
├── alerax.log                        # Run metadata and convergence log
├── ccpdim.txt                        # CCP dimensions per family
├── fractionMissing.txt               # Per-species missing data fractions
├── perSpeciesCoverage.txt            # Per-species family coverage
├── per_fam_likelihoods.txt           # Log-likelihood per family
├── model_parameters/
│   └── model_parameters.txt          # Optimised DTL rates per gene
├── species_trees/
│   ├── inferred_species_tree.newick  # Final optimised species tree
│   └── starting_species_tree.newick  # Input species tree
└── reconciliations/
    ├── transfers.txt                 # Global transfer events
    ├── perspecies_eventcount.txt      # Global per-species event counts
    ├── origins/                       # Origination probability per node
    ├── all/                           # Per-family, per-sample results
    │   ├── {FAMILY}.newick           # All sampled gene trees
    │   ├── {FAMILY}.rec_uml          # Annotated reconciled tree
    │   ├── {FAMILY}_{N}.xml          # Detailed XML reconciliation (N=0..99)
    │   ├── {FAMILY}_eventcount_{N}.txt
    │   ├── {FAMILY}_transfers_{N}.txt
    │   └── {FAMILY}_perspecies_eventcount_{N}.txt
    └── summaries/                     # Per-family summaries across samples
        ├── {FAMILY}_consensus_50.newick
        ├── {FAMILY}_transfers.txt
        └── {FAMILY}_perspecies_eventcount.txt
```

#### AleRax run-level analysis

```python
from reconparser import AleRaxRun

run = AleRaxRun("output_dir")

# All available families
families = run.get_family_names()
# ['K00192', 'K00193', ...]

# Per-species events aggregated across all families
perspecies = run.get_perspecies_events()

# Origination probabilities for a species-tree node
origin = run.get_origin("Node_a1001_a1048_0")
# {'vertical': 0.1, 'b93910': 0.01, ...}

# Coverage and missing data
coverage = run.get_per_species_coverage()
missing = run.get_fraction_missing()
```

#### AleRax per-family analysis

```python
fam = run.get_family("K00192")

# How many sampled gene trees?
n = fam.get_num_samples()  # e.g. 100

# Per-species event counts (averaged across samples)
perspecies = fam.get_perspecies_events()

# Transfer events for a specific sample
sample_transfers = fam.get_sample_transfers(sample=0)

# Per-species events for a specific sample
sample_events = fam.get_sample_perspecies_events(sample=0)

# Check which files exist for this family
fam.files_exist()
# {'newick': True, 'rec_uml': True, 'consensus_tree': True, ...}
```

## API Reference

### ALEParser

```python
ALEParser(base_path: str | Path)
```

Parser for classic ALE output. Accepts a base path (e.g., `"results.ale"`) or a specific file path (e.g., `"results.ale.ucons_tree"`) — the base path is extracted automatically.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_consensus_tree()` | `ete3.Tree` | Consensus gene tree from `.ucons_tree` |
| `get_transfers()` | `pd.DataFrame` | Transfer events (columns: `from`, `to`, `freq`) |
| `get_transfers_as_dict_list()` | `List[Dict]` | Transfers as list of dicts |
| `get_reconciled_tree()` | `ete3.Tree` | Species tree from `.uml_rec` |
| `get_reconciled_gene_trees()` | `List[ete3.Tree]` | All reconciled gene trees (with D@/T@ annotations) |
| `get_ml_rates()` | `Dict[str, float]` | ML rates: `duplications`, `transfers`, `losses` |
| `get_log_likelihood()` | `float` | Log-likelihood of reconciliation |
| `get_summary_statistics()` | `Dict[str, float]` | Total event counts across all branches |
| `get_branch_statistics()` | `pd.DataFrame` | Per-branch event statistics |
| `get_all_data()` | `Dict` | Everything in one call (gracefully handles missing files) |
| `files_exist()` | `Dict[str, bool]` | Check which output files are present |

### AleRaxRun

```python
AleRaxRun(output_dir: str | Path)
```

Run-level parser for an AleRax output directory. All data is loaded lazily.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_run_info()` | `Dict` | Run metadata from `alerax.log` (version, model, counts, etc.) |
| `get_family_names()` | `List[str]` | Sorted list of gene family identifiers |
| `get_family(name)` | `AleRaxFamily` | Lazy per-family parser |
| `get_families()` | `Dict[str, AleRaxFamily]` | All family parsers keyed by name |
| `get_species_tree()` | `ete3.Tree` | Inferred (optimised) species tree |
| `get_starting_species_tree()` | `ete3.Tree` | Input species tree |
| `get_model_parameters()` | `pd.DataFrame` | DTL rates per gene (`gene`, `dup_rate`, `loss_rate`, `transfer_rate`) |
| `get_per_family_likelihoods()` | `pd.DataFrame` | Per-family log-likelihoods (`family`, `log_likelihood`) |
| `get_total_log_likelihood()` | `float` | Sum of per-family log-likelihoods |
| `get_transfers()` | `pd.DataFrame` | Global transfers (`from`, `to`, `score`) |
| `get_perspecies_events()` | `pd.DataFrame` | Global per-species event counts |
| `get_origin(node_name)` | `Dict[str, float]` | Origination probabilities for a species-tree node |
| `get_origin_node_names()` | `List[str]` | Node names with origin files |
| `get_fraction_missing()` | `pd.DataFrame` | Per-species missing data fraction |
| `get_per_species_coverage()` | `pd.DataFrame` | Per-species family coverage |
| `get_ccp_dimensions()` | `pd.DataFrame` | CCP dimensions per family |
| `files_exist()` | `Dict[str, bool]` | Check which run-level files/directories exist |

### AleRaxFamily

```python
AleRaxFamily(family_name: str, output_dir: str | Path)
```

Per-family parser for AleRax output. Usually obtained via `AleRaxRun.get_family()`. All data is loaded lazily and cached.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_num_samples()` | `int` | Number of sampled gene trees |
| `get_consensus_tree()` | `ete3.Tree` | Majority-rule consensus gene tree (with support values) |
| `get_sampled_gene_tree(i)` | `ete3.Tree` | Single sampled gene tree (memory-efficient) |
| `get_sampled_gene_trees()` | `List[ete3.Tree]` | All sampled gene trees (loads all into memory) |
| `get_reconciled_gene_tree()` | `ete3.Tree` | Annotated reconciled tree (`.rec_uml`) with T@/D@/S labels |
| `get_transfers()` | `pd.DataFrame` | Averaged transfers (`from`, `to`, `freq`) |
| `get_transfers_as_dict_list()` | `List[Dict]` | Averaged transfers as list of dicts |
| `get_perspecies_events()` | `pd.DataFrame` | Per-species events averaged across samples |
| `get_event_counts(sample)` | `Dict[str, int]` | Event counts for one sample (S, SL, D, DL, T, TL, L, Leaf) |
| `get_all_event_counts()` | `pd.DataFrame` | Event counts for all samples (one row per sample) |
| `get_sample_transfers(sample)` | `pd.DataFrame` | Transfer events for one sample |
| `get_sample_perspecies_events(sample)` | `pd.DataFrame` | Per-species events for one sample |
| `files_exist()` | `Dict[str, bool]` | Check which output files exist for this family |

## Design Principles

All parsers in reconparser follow the same conventions:

- **Lazy loading** — Data is parsed only when you call the corresponding getter. This is especially important for AleRax output, which can contain millions of files across thousands of gene families.
- **Caching** — Once parsed, data is stored in memory so repeated calls are instant.
- **pandas DataFrames** — Tabular data is returned as DataFrames for easy filtering, grouping, and export.
- **ete3 Trees** — Phylogenetic trees are returned as `ete3.Tree` objects for traversal, annotation, and visualisation.
- **FileNotFoundError** — Missing files raise clear errors with the expected path.

## File Format Support

| Software | Version | Status |
|----------|---------|--------|
| ALE | v0.4, v1.0 | Fully supported |
| AleRax | v1.2+ | Fully supported (except XML reconciliations) |

Future releases will include parsers for additional reconciliation tools.

## Development

```bash
cd reconparser
pip install -e ".[dev]"
pytest tests/
```

### Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/new-parser`)
3. Commit your changes (`git commit -am 'Add parser for ToolX'`)
4. Push to the branch (`git push origin feature/new-parser`)
5. Open a Pull Request

## Citation

If you use reconparser in your research, please cite:

```
@software{reconparser2024,
  author = {Davin, Adrian A.},
  title = {reconparser: A Python library for parsing phylogenetic reconciliation outputs},
  year = {2024},
  url = {https://github.com/aadavin/reconparser}
}
```

And please also cite the original reconciliation software you used (e.g., ALE, AleRax).

## License

MIT License — see [LICENSE](LICENSE) file for details.

## Links (not available yet)

- **GitHub**: https://github.com/aadavin/reconparser
- **PyPI**: https://pypi.org/project/reconparser/
- **Issues**: https://github.com/aadavin/reconparser/issues

## Acknowledgments

reconparser is developed and maintained by Adrian A. Davin.

Special thanks to the developers of:
- **ALE** — Szöllősi GJ, et al. (2013)
- **AleRax** — Morel B, et al. (2022)
- **ete3** — Huerta-Cepas J, et al. (2016)
