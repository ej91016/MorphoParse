# MorphoParse

**A Python package and command-line tool for parsing, partitioning, and weighting morphological character matrices for phylogenetic analysis**

![GitHub release](https://img.shields.io/github/tag/ej91016/MorphoParse.svg)
<p align="left">
  <img src="mp_logo.png" alt="MorphoParse Logo" width="300"/>
</p>

## 🦣 Who Should Use MorphoParse
MorphoParse is for you if:
- You're trying to input a **morphological matrix** into phylogenetic software and:
  - You know what’s wrong, but want a tool to fix it automatically, or
  - You don’t know what’s wrong and need help identifying and resolving issues
- You need a tool in your analysis pipeline to:
  - Standardize **matrix format**
  - Standardize **polymorphic encodings**
  - Convert between matrix formats
- You want to **model with correct state space** in phylogenetic inference

## ✨ Features
- Supports **FASTA**, **PHYLIP**, **NEXUS**, and **TNT** formats
- Normalizes taxon names by stripping quotes and replacing spaces/symbols with `_`
- Modify data for better software compatibility
  - Optional removal of **polymorphic encodings** (e.g., `{01}`, `[0 1]`) and replaces them with `?`
    - These will errors in **RAxML** and are ignored in **IQ-TREE**
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


## 🧱 Requirements

- Python **3.7+**
- [Biopython](https://biopython.org/)

Install via pip:

```bash
pip install biopython
```

## 📦 Installation

Install the latest **stable** or **pre-release version** (recommended):

```bash
pip install git+https://github.com/ej91016/MorphoParse.git@v0.1.0b0
```

Install the latest **development version** directly from GitHub:

```bash
pip install git+https://github.com/ej91016/MorphoParse.git
```



## ⚙️ Command-Line Options

| Flag                 | Description                                                           |
| -------------------- | -------------------------------------------------------------------   |
| `-i`, `--input`      | Input morphological matrix (required)                                 |
| `-o`, `--output`     | Output prefix (default: input filename without extension)             |
| `-f`, `--format`     | Input format: `fasta`,`phylip`,`nexus`,`tnt` (default: auto-detection)|
| `-g`, `--out_format` | Output format (default: same as input)                                |
| `-p`, `--poly`       | Keep polymorphic encodings (avoid if using maximum likelihood)        |
| `-r`, `--remap`      | Enable all of: remove missing, remove mono, reorder states            |
| `--remove-missing`   | Remove characters with only `?` or `-`                                |
| `--remove-mono`      | Remove characters with only one unambiguous state                     |
| `--reorder`          | Renumber states to 0,1,... per site                                   |
| `-a`, `--asc`        | Apply ASC correction to partition models                              |
| `-s`, `--software`   | Software: `raxml`,`raxmlng`,`iqtree` (default: `raxml`)               |
| `--paup`             | Generate SSA weights for PAUP\*                                       |
| `--tnt`              | Generate SSA weights for TNT                                          |
| `--version`          | Show version info and exit                                            |

---
## 📊 MorphoParse Flag Recommendations by Analysis Type

- **SSA**: State-space-aware (the "correct" setup when doing just one analysis)
- **SSM**: State-space-misspecified (a special case of False-space)
- **False-space**: Purposely misspecify the state space


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
morphoparse -i example.nexus -g phylip -r -a --paup
```

This will:

- Parse `example.nexus` as a **NEXUS** file
- Replace polymorphic encodings with `?`
- Produce output files in the current directory:

  - `example_mparse.phy`         – Parsed and polished data
  - `example_remap.txt`          – Remapping details
  - `example_raxml.models`       – ASC-corrected partitions for RAxML
  - `example_paup_weights.txt`   – PAUP\* weights
<br><br>

Example data can be found in [`examples/`](https://github.com/ej91016/MorphoParse/tree/main/examples)


## 📝 Troubleshooting

To ensure proper SSA setup for IQ-TREE:
- Choose `PHYLIP` or `FASTA` as output format
- Use only `-p *_iqtree.nex` in IQ-TREE (not `-s`)
  - Defines partitions & linked matrix files
  - Assume IQ-TREE will be called from the same directory as MorphoParse — edit paths if needed
- If you are running IQ-TREE2, replace MK with JC2 for the binary partition
- If you see the error `ERROR: Cannot concatenate sub-alignments of different #states`, add the flag `-keep-ident`
- If IQ-TREE end prematurely after showing initial log-Likelihood:
  - IQ-TREE appears to have trouble with `MK+ASC` (with or without `+G`)
  - *Fix (for now)*: Remove `+ASC` from the models (or avoid `-a` in MorphoParse)
  - Alternatively, consider using RAxML or RAxML-NG

---

## 🪪 License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)
