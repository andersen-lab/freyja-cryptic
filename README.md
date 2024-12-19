# Freyja-Cryptic

Detect cryptic variant circulation in wastewater sequencing data. Beginning with ![Freyja covariants](https://andersen-lab.github.io/Freyja/src/usage/covariants.html) outputs, this tool will query each mutation cluster occurring on the same read pair to determine if the cluster has been previously detected in clinical sequencing.

## Installation

```bash
git clone https://github.com/andersen-lab/freyja-cryptic.git && cd freyja-cryptic

mamba create -n freyja-cryptic
mamba env update -n freyja-crypic -f environment.yml
mamba activate freyja-cryptic
```

## Usage

```bash
python freyja-cryptic.py --covariants_dir /path/to/covariants --output_dir /path/to/output --metadata /path/to/metadata.csv --max_clinical_count (default=10)
```

Upon running the script, a prompt will appear asking to authenticate with GISAID. This is needed to access the clinical data.