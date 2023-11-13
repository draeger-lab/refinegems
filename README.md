[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fdraeger-lab%2Frefinegems%2Fmain%2Fpyproject.toml)
[![Documentation Status](https://readthedocs.org/projects/refinegems/badge/?version=latest)](https://refinegems.readthedocs.io/en/latest/?badge=latest)
![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/refinegems?logo=github&label=refineGEMs&color=B4A069&style=flat-square)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/draeger-lab/refinegems/main)
![Repo Size](https://img.shields.io/github/repo-size/draeger-lab/refinegems)
![GitHub all releases](https://img.shields.io/github/downloads/draeger-lab/refinegems/total?logo=github&label=GitHub%20downloads)
[![PyPI version](https://img.shields.io/pypi/v/refinegems?logo=pypi&label=PyPI%20package&color=neongreen)](https://pypi.org/project/refineGEMs/)
![PyPI - Format](https://img.shields.io/pypi/format/refinegems?)
[![PyPI downloads](https://img.shields.io/pypi/dm/refinegems.svg?logo=pypi&label=PyPI%20downloads)](https://pypistats.org/packages/refinegems)  
[![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8270303-B4A069?style=flat-square&logo=zenodo&logoColor=white)](https://zenodo.org/badge/latestdoi/359867657)  
[![Frontiers DOI](https://img.shields.io/badge/Frontiers%20DOI-10.3389%2Ffbinf.2023.1214074-B4A069?style=flat-square)](https://www.frontiersin.org/articles/10.3389/fbinf.2023.1214074/full)

<p align="center">
<img src="https://github.com/draeger-lab/refinegems/raw/main/docs/source/images/refineGEMs_logo.png" height="200"/>
</p>

# refineGEMs
`refineGEMs` is a python package intended to help with the curation of genome-scale metabolic models (GEMS). </br>
The documentation can be found [here](https://refinegems.readthedocs.io/en/latest/).

## Table of contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [How to cite](#how-to-cite)
4. [Repositories using refineGEMs](#repositories-using-refinegems)

## Overview

Currently `refineGEMs` can be used for the investigation of a GEM, it can complete the following tasks:

- loading GEMs with `COBRApy` and `libSBML`
- report number of metabolites, reactions and genes
- report orphaned, deadends and disconnected metabolites
- report mass and charge unbalanced reactions
- report [Memote](https://memote.readthedocs.io/en/latest/index.html) score
- compare the genes present in the model to the genes found in:
  - the [KEGG](https://www.genome.jp/kegg/kegg1.html) Database (Note: This requires the GFF file and the KEGG identifier of your organism.)
  - Or the [BioCyc](https://biocyc.org) Database (Note: This requires that a database entry for your organism exists in BioCyc.)
- compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the [ModelSEED](https://modelseed.org/) Database.

Other applications of `refineGEMs` to curate a given model include: 

- The correction of a model created with [CarveMe](https://github.com/cdanielmachado/carveme) v1.5.1 or v1.5.2 (for example moving all relevant information from the notes to the annotation field or automatically annotating the GeneProduct section of the model with the respective NCBI gene/protein identifiers from the GeneProduct identifiers),
- The addition of [KEGG](https://www.genome.jp/kegg/kegg1.html) Pathways as Groups (using the [libSBML](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html) Groups Plugin),
- Updating the SBO-Term annotations based on [SBOannotator](https://github.com/draeger-lab/SBOannotator),
- Updating the annotation of metabolites and extending the model with reactions (for the purpose of filling gaps) based on a table filled by the user `data/manual_annotations.xlsx` (Note: This only works when the structure of the [example Excel file](https://github.com/draeger-lab/refinegems/blob/5eac900d9848b5ae5faf0055db72a986e7ba64e8/data/manual_curation.xlsx) is used.),
- And extending the model with all information surrounding reactions including the corresponding GeneProducts and metabolites by filling in the table `data/modelName_gapfill_analysis_date_example.xlsx` (Note: This also only works when the structure of the [example Excel file](https://github.com/draeger-lab/refinegems/blob/5eac900d9848b5ae5faf0055db72a986e7ba64e8/data/modelName_gapfill_analysis_date_example.xlsx) is used).

## Installation

You can install `refineGEMs` via pip:

```bash
pip install refineGEMs

```

or to a local conda environment where `refineGEMs` is distributed via this GitHub repository and all dependencies are denoted in the `setup.py` file:

```bash
# clone or pull the latest source code
git clone https://github.com/draeger-lab/refinegems.git
cd refinegems

conda create -n <EnvName> python=3.9

conda activate <EnvName>

# check that pip comes from <EnvName>
which pip

pip install .

```

## How to cite
When using `refineGEMs`, please cite the latest publication:

Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas Dräger. 
Genome-scale metabolic models consistently predict in vitro characteristics of Corynebacterium
striatum. Front. Bioinform., oct 2023. [doi:10.3389/fbinf.2023.1214074](https://doi.org/10.3389/fbinf.2023.1214074).

## Repositories using refineGEMs
- [C_striatum_GEMs](https://github.com/draeger-lab/C_striatum_GEMs)
- draeger-lab/Shaemolyticus - `private`
- draeger-lab/Ssanguinis - `private`
