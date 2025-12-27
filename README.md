[![license](https://img.shields.io/badge/license-GPL%20v3-black.svg?style=flat-square)](LICENSE)

# selscape-analysis

## Introduction

This repository contains a Snakemake workflow designed to reproduce the results of analyses for natural selection using [selscape](https://github.com/xin-huang/selscape). The workflow has been tested on Rocky Linux 9.6 (Blue Onyx) using the [Life Science Compute Cluster](https://lisc.univie.ac.at/) at the University of Vienna and the [Multi-Site Computer Austria](https://docs.asc.ac.at/systems/musica.html) of Austrian Scientific Computing.

## Usage

1. Install [miniforge](https://github.com/conda-forge/miniforge/releases). [Miniforge3 (version 25.9.1)](https://github.com/conda-forge/miniforge/releases/download/25.9.1-0/Miniforge3-25.9.1-0-Linux-x86_64.sh) was used for analysis.

2. Clone this repository:

```
git clone https://github.com/xin-huang/selscape-analysis
cd selscape-analysis
```

3. Create the environment:

```
mamba env create -f env.yaml
```

4. Activate the environment:

```
mamba activate selscape-analysis
```

5. Download data and `selscape`:

```
snakemake -c 1 -s download_resources.smk
```

6. Run the analysis locally:

```
snakemake -c 1 --configfile config/main.yaml
```

7. Run the analysis on HPC:

```
snakemake -c 1 --configfile config/main.yaml --profile config/slurm/lisc
```

Users should adjust the resource parameters in each Snakemake file to match their cluster settings and modify the `config.yaml` file in `config/slurm` according to their job scheduler.
Before analysis, users should manually download [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and place it in `resources/tools`.
