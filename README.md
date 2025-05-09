<p align="center">
<img src="https://github.com/draeger-lab/refinegems/raw/main/docs/source/images/refineGEMs_logo.png" height="200"/>
</p>

| Topic | Badge(s) |
| :--- | :---- |
| General | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fdraeger-lab%2Frefinegems%2Fmain%2Fpyproject.toml) [![Documentation Status](https://readthedocs.org/projects/refinegems/badge/?version=latest)](https://refinegems.readthedocs.io/en/latest/?badge=latest) ![Repo Size](https://img.shields.io/github/repo-size/draeger-lab/refinegems) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/draeger-lab/refinegems/main) |
| GitHub release | ![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/refinegems?logo=github&label=refineGEMs&color=B4A069&style=flat-square) ![GitHub all releases](https://img.shields.io/github/downloads/draeger-lab/refinegems/total?logo=github&label=GitHub%20downloads) |
| PyPI | [![PyPI version](https://img.shields.io/pypi/v/refinegems?logo=pypi&label=PyPI%20package&color=neongreen)](https://pypi.org/project/refineGEMs/) [![PyPI downloads](https://img.shields.io/pypi/dm/refinegems.svg?logo=pypi&label=PyPI%20downloads)](https://pypistats.org/packages/refinegems) ![PyPI - Format](https://img.shields.io/pypi/format/refinegems?) |
| Compliance | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=32102&r=133) [![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10532/badge)](https://www.bestpractices.dev/projects/10532) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| References | [![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8270303-B4A069?style=flat-square&logo=zenodo&logoColor=white)](https://zenodo.org/badge/latestdoi/359867657) [![Frontiers DOI](https://img.shields.io/badge/Frontiers%20DOI-10.3389%2Ffbinf.2023.1214074-B4A069?style=flat-square)](https://www.frontiersin.org/articles/10.3389/fbinf.2023.1214074/full) |


# refineGEMs
`refineGEMs` is a python package intended to help with the curation of genome-scale metabolic models (GEMS). </br>

> [!WARNING]
> 🚧 The documentation is currently under heavy-rework!
> The documentation can be found [here](https://refinegems.readthedocs.io/en/latest/).🚧

## Table of contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [How to cite](#how-to-cite)
4. [Repositories using refineGEMs](#repositories-using-refinegems)

## Overview

Currently `refineGEMs` can be used for the investigation of a genome-scale metabolic model (GEM)/multiple GEMs, it can complete the following tasks:

- Loading GEMs with `COBRApy` and `libSBML`
- Report and visualise number of metabolites, reactions and genes
- Report orphaned, deadends and disconnected metabolites
- Report mass and charge unbalanced reactions
- Report the [Memote](https://memote.readthedocs.io/en/latest/index.html) score and provide a whole MEMOTE report
- Find and fill gaps automatically via databases like KEGG, BioCyc, SwissProt or a user-defined database
- Compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the [ModelSEED](https://modelseed.org/) Database.

Other applications of `refineGEMs` to curate a given model include: 

- The correction of a model created with [CarveMe](https://github.com/cdanielmachado/carveme) v1.5.1 or v1.5.2 (for example moving all relevant information from the notes to the annotation field or automatically annotating the GeneProduct section of the model with the respective NCBI gene/protein identifiers from the GeneProduct identifiers)
- The addition of [KEGG](https://www.genome.jp/kegg/kegg1.html) Pathways as Groups (using the [libSBML](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html) Groups Plugin)
- Updating the SBO-Term annotations based on [SBOannotator](https://github.com/draeger-lab/SBOannotator)


## Installation

You can install `refineGEMs` via pip:

```bash
pip install refineGEMs

```

or to a local conda environment where `refineGEMs` is distributed via this GitHub repository and all dependencies are denoted in the `pyproject.toml` file:

```bash
# clone or pull the latest source code
git clone https://github.com/draeger-lab/refinegems.git
cd refinegems

conda create -n <EnvName> python=3.10 (at least but < 3.13)

conda activate <EnvName>

# check that pip comes from <EnvName>
which pip

pip install .

```

> [!CAUTION]
> ``refineGEMs`` depends on the tools [MCC](https://github.com/Biomathsys/MassChargeCuration) and 
> [BOFdat](https://github.com/draeger-lab/BOFdat) which cannot directly be installed via PyPI or the `pyproject.toml`. 
> Please install both tools before using ``refineGEMs``:
>
> ```bash
> # For MCC, until hot fix is merged into main:
> pip install "masschargecuration@git+https://github.com/Biomathsys/MassChargeCuration"
>
> # For BOFdat, our fork with hot fix(es):
> pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"
>
> ```

## How to cite
When using `refineGEMs`, please cite the latest publication:

Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas Dräger. 
Genome-scale metabolic models consistently predict in vitro characteristics of Corynebacterium
striatum. Front. Bioinform., oct 2023. [doi:10.3389/fbinf.2023.1214074](https://doi.org/10.3389/fbinf.2023.1214074).

## Repositories using refineGEMs
- [C_striatum_GEMs](https://github.com/draeger-lab/C_striatum_GEMs)
- draeger-lab/Shaemolyticus - `private`
- draeger-lab/Ssanguinis - `private`
