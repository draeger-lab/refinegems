[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7915766.svg)](https://doi.org/10.5281/zenodo.7915766) [![Documentation Status](https://readthedocs.org/projects/refinegems/badge/?version=latest)](https://refinegems.readthedocs.io/en/latest/?badge=latest) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
<img src="https://github.com/draeger-lab/refinegems/raw/main/docs/source/images/refineGEMs_logo.png" height="200"/>
</p>

# refineGEMs
`refineGEMs` is a python package inteded to help with the curation of genome-scale metabolic models (GEMS).

## Documentation
The docs can be found [here](https://refinegems.readthedocs.io/en/latest/).

## Overview

Currently `refineGEMs` can be used for the investigation of a GEM, it can complete the following tasks:
- loading GEMS with `cobrapy` and `libSBML`
- report number of metabolites, reactions and genes
- report orphaned, deadends and disconnected metabolites
- report mass and charge unbalanced reactions
- report [Memote](https://memote.readthedocs.io/en/latest/index.html) score
- compare the genes present in the model to the genes found in the [KEGG](https://www.genome.jp/kegg/kegg1.html) Database (Note: this requires a gff file of your organism and the KEGG identifier of your organism)
- compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the [ModelSEED](https://modelseed.org/) Database

Other applications of `refineGEMs` include curation of a given model these include:
- correction of a model created with [CarveMe](https://github.com/cdanielmachado/carveme) v.1.5.1 (for example moving all relevant information from the notes to the annotation field) this includes automated annotation of NCBI genes to the GeneProtein section of the model
- addition of [KEGG](https://www.genome.jp/kegg/kegg1.html) Pathways as Groups (using the [libSBML](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html) Groups Plugin)
- SBO-Term annotation based on a script by Elisabeth Fritze
- annotation of metabolites based using a table created by the user `data/manual_annotations.xlsx`

## Installation

You can install `refineGEMs` via pip:

```bash
pip install refinegems

```

or to a local conda environment where `refineGEMs` is distributed via this github repository and all dependencies are denoted in the `setup.py` file:

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
