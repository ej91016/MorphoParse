# MorphoParse

**A Python package and command-line tool for parsing, partitioning, and weighting morphological character matrices for phylogenetic analysis**

## ✨ Features

- Supports **FASTA**, **PHYLIP**, **NEXUS**, and **TNT** input formats
- Removes **polymorphic encodings** (e.g., `{01}`, `[0,1]`, `(0/1)`) and replaces them with `?`

  - These cause errors in **RAxML** and are ignored in **IQ-TREE**
- Normalizes taxon names by stripping quotes and replacing spaces/symbols with `_`
- Optional **state remapping** (e.g., `0,1,4,5` → `0,1,2,3`)

  - Required by some tools (e.g., **RAxML**)
- Optional removal of uninformative or missing-only characters
- Outputs clean alignment in selected format
- Generates **partition model files** compatible with:
  - RAxML
  - RAxML-NG
  - IQ-TREE
- Generates **state-space-aware weighting (SSA)**:
  - `SSA weight = ln(r)`, where `r` = number of observed states
  - Supports **PAUP\*** and **TNT**

    - Weights are implemented as `round(10 × ln(r))`
    - Both programs accept only integer weights
  - See [Huang (2025, preprint)](https://doi.org/10.1101/2025.04.22.650124) for rationale

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

| Flag                 | Description                                                         |
| -------------------- | ------------------------------------------------------------------- |
| `-i`, `--input`      | Input morphological matrix (required)                               |
| `-o`, `--output`     | Output prefix (default: input filename without extension)           |
| `-f`, `--format`     | Input format: `fasta`, `phylip`, `nexus`, `tnt` (default: `phylip`) |
| `-g`, `--out_format` | Output format (default: same as input)                              |
| `-r`, `--remap`      | Enable all of: remove missing, remove mono, reorder states          |
| `--remove-missing`   | Remove characters with only `?` or `-`                              |
| `--remove-mono`      | Remove characters with only one unambiguous state                   |
| `--reorder`          | Renumber states to 0,1,... per site                                 |
| `-a`, `--asc`        | Apply ASC correction to partition models                            |
| `-s`, `--software`   | Software: `raxml`,`raxmlng`,`iqtree` (default: `raxml`)             |
| `-p`, `--paup`       | Generate SSA weights for PAUP\*                                     |
| `-t`, `--tnt`        | Generate SSA weights for TNT                                        |
| `--version`          | Show version info and exit                                          |

---
## 📊 MorphoParse Flag Recommendations by Analysis Type

- **SSA**: State-space-aware (the "correct" setup when doing just one analysis)
- **SSM**: State-space-misspecified (a special case of False-space)
- **False-space**: Purposely misspecify the state space

Please refer to [Huang (2025, preprint)](https://doi.org/10.1101/2025.04.22.650124) for more details

**Legend:** ✅ Required  ⭐ Recommended  🔘 Optional  ✘ No

| Analysis Type        | Remap               | Remove Mono  | Remove Missing  | ASC | Partition File  | Notes                      | Tools    |
| -------------------- | ------------------- | ------------ | --------------- | --- | --------------- | -------------------------- | -------- |
| SSA (with invariant) | ⭐                  | ✘           | ✅              | ✘   | ✅             | Recommended model          | All      |
| SSA                  | ⭐                  | ✅          | ✅              | ✅  | ✅             | No invariant sites allowed | All      |
| SSM (with invariant) | ⭐                  | ✘           | ✅              | ✘   | ✘              | Default model              | All      |
| SSM                  | ⭐                  | ✅          | ✅              | ✅  | ✘              | No invariant sites allowed | All      |
| FS (padding)         | ✘ (destroy padding) | 🔘          | ✅              | ✘   | ✘              | Incompatible with ASC      | All      |
| FS (override)        | ⭐                  | 🔘          | ✅              | 🔘  | 🔘             | Specify with MULTI`x`\_MK  | RAxML-NG |

---

## 🥪 Usage

Example usage:
```bash
morphoparse -i example.nexus -f nexus -g phylip -r -a -p
```

This will:

- Parse the **NEXUS** file: `example.nexus`
- Produce output files in the current directory:

  - `example_clean.phy`          – Cleaned and remapped data
  - `example_remap.txt`          – Remapping details
  - `example_raxml.models`       – ASC-corrected partitions
  - `example_paup_weights.txt`   – PAUP\* weights
<br><br>

Example data can be found in [`examples/`](https://github.com/ej91016/MorphoParse/tree/main/examples)
<br><br>

📝 Note: Proper SSA setup for IQ-TREE 
- Choose `PHYLIP` or `NEXUS` as output format
- Use only `-p *_iqtree.nex` in IQ-TREE (not `-s`)
  - Defines partitions & linked matrix files
  - Assume IQ-TREE will be called from the same directory as MorphoParse — edit paths if needed
- If IQ-TREE end prematurely after showing initial log-Likelihood:
  - IQ-TREE appears to have trouble with `MK+ASC` (with or without `+G`)
  - *Fix (for now)*: Remove `+ASC` from the models (or avoid `-a` in MorphoParse)
  - Alternatively, consider using RAxML or RAxML-NG

---

## 📖 Citation

If you use **MorphoParse**, please cite:

> Huang EJ (2025). *State Space Misspecification in Morphological Phylogenetics: A Pitfall for Models and Parsimony Alike*.
> bioRxiv. [https://doi.org/10.1101/2025.04.22.650124](https://doi.org/10.1101/2025.04.22.650124)

## 🪪 License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)
