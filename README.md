# MorphoParse

**A Python package and command-line tool for parsing, partitioning, and weighting morphological character matrices for phylogenetic analysis**

## ✨ Features

- Supports **FASTA**, **PHYLIP**, **NEXUS**, and **TNT** input formats
- Removes **polymorphic encodings** (e.g., `{01}`, `[0,1]`, `(0/1)`) and replaces them with `?`
  - These cause errors in **RAxML** and are ignored in **IQ-TREE 2**
- Normalizes taxon names by stripping quotes and replacing spaces/symbols with `_`
- Optional **state remapping** (e.g., `0,1,4,5` → `0,1,2,3`)
  - Required by some tools (e.g., **RAxML**)
- Optional removal of uninformative or missing-only characters
- Outputs clean alignment in selected format
- Generates **partition model files** compatible with:
  - RAxML
  - RAxML-NG
  - IQ-TREE
- Supports **state-space-aware weighting (SSA)** for:
  - **PAUP\***
  - **TNT**
  - Uses `SSA weight = round(10 × ln(r))`, where `r = number of observed states`

## 🧱 Requirements

- Python **3.7+**
- [Biopython](https://biopython.org/)

Install via pip:

```bash
pip install biopython
```

## 📦 Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/ej91016/MorphoParse.git
```

## ⚙️ Command-Line Options

| Flag | Description |
|------|-------------|
| `-i`, `--input` | Input morphological matrix (required) |
| `-o`, `--output` | Output prefix (default: input filename without extension) |
| `-f`, `--format` | Input format: `fasta`, `phylip`, `nexus`, `tnt` (default: `phylip`) |
| `-g`, `--out_format` | Output format (default: same as input) |
| `-r`, `--remap` | Enable all of: remove missing, remove mono, reorder states |
| `--remove-missing` | Remove characters with only `?` or `-` |
| `--remove-mono` | Remove characters with only one unambiguous state |
| `--reorder` | Renumber states to 0,1,... per site |
| `-a`, `--asc` | Apply ASC correction to partition models |
| `-n`, `--raxmlng` | Use RAxML-NG partition format (default: RAxML) |
| `-p`, `--paup` | Generate SSA weights for PAUP\* |
| `-t`, `--tnt` | Generate SSA weights for TNT |
| `--version` | Show version info and exit |

## 🧪 Example Usage

```bash
morphoparse -i example.nexus -f nexus -g phylip -r -a -p
```

This will:
- Parse the **NEXUS** file: `example.nexus`
- Produce output files in the current directory:
  - `example_clean.phy`          – Cleaned and remapped data
  - `example_remap.txt`          – Remapping details
  - `example_raxml.models`       – ASC-corrected partitions
  - `example_paup_weights.txt`   – PAUP* weights

> **📝 Note:** Example data can be found in the [`examples/`](https://github.com/ej91016/MorphoParse/tree/main/examples) directory of the [GitHub repository](https://github.com/ej91016/MorphoParse).

---

## 📊 MorphoParse Flag Recommendations by Analysis Type

(paper describing analysis type coming soon...)

**Legend:** ✅ Required  ⭐ Recommended  🔘 Optional  ✘ No

| Analysis Type | Remap | Remove Mono | Remove Missing | ASC | Partition File | Notes | Tools |
|-------------------|--------|--------------|----------------|-----|----------------|-------|-------|
| SSM (invariant) | ⭐ | ✘ | ✅ | ✘ | ✘ | Default model | All |
| SSM  | ⭐ | ✅ | ✅ | ✅ | ✘ | No invariant sites allowed | All |
| SSA (invariant) | ⭐ | ✘ | ✅ | ✘ | ✅ | IQ-TREE accepts RAxML-style | All |
| SSA | ⭐ | ✅ | ✅ | ✅ | ✅ | No invariant sites allowed | All |
| FS (padding) | ✘ (destroy padding) | 🔘 | ✅ | ✘ | ✘ | Incompatible with ASC | All |
| FS (override) | ⭐ | 🔘 | ✅ | 🔘 | 🔘 | Specify with MULTI`x`_MK| RAxML-NG |

---

## 📖 Citation

If you use **MorphoParse**, please cite:

> EJ Huang (2025). State Space Misspecification in Morphological Phylogenetics: A Pitfall for Models and Parsimony Alike. bioRxiv. https://doi.org/10.1101/2025.04.22.650124

## 🪪 License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).
