<p align="center">
<img src="https://github.com/draeger-lab/refinegems/raw/main/docs/source/images/logos/full-nb.png" height="200"/>
</p>

| Topic | Badge(s) |
| :--- | :---- |
| General | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fdraeger-lab%2Frefinegems%2Fmain%2Fpyproject.toml) [![Documentation Status](https://readthedocs.org/projects/refinegems/badge/?version=latest)](https://refinegems.readthedocs.io/en/latest/?badge=latest) ![Repo Size](https://img.shields.io/github/repo-size/draeger-lab/refinegems) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/draeger-lab/refinegems/main) |
| GitHub release | ![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/refinegems?logo=github&label=refineGEMs&color=B4A069&style=flat-square) ![GitHub all releases](https://img.shields.io/github/downloads/draeger-lab/refinegems/total?logo=github&label=GitHub%20downloads) |
| PyPI | [![PyPI version](https://img.shields.io/pypi/v/refinegems?logo=pypi&label=PyPI%20package&color=neongreen)](https://pypi.org/project/refineGEMs/) ![PyPI - Format](https://img.shields.io/pypi/format/refinegems?) |
| Compliance | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=32102&r=133) [![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10532/badge)](https://www.bestpractices.dev/projects/10532) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| References | [![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8270303-B4A069?style=flat-square&logo=zenodo&logoColor=white)](https://zenodo.org/badge/latestdoi/359867657) [![Frontiers DOI](https://img.shields.io/badge/Frontiers%20DOI-10.3389%2Ffbinf.2023.1214074-B4A069?style=flat-square)](https://www.frontiersin.org/articles/10.3389/fbinf.2023.1214074/full) |


# refineGEMs
`refineGEMs` is a python package intended to help with the curation of genome-scale metabolic models (GEMS). </br>

<!-- TOC -->

- [Overview](#overview)
- [Installation](#installation)
    - [Via pip](#pypi-via-pip)
    - [Via Docker](#docker-via-docker)
- [License](#license)
- [How to cite](#how-to-cite)
- [Repositories using refineGEMs](#repositories-using-refinegems)

<!-- /TOC -->

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
The toolbox ``refineGEMs``can be installed via pip or via Docker.

### ![pypi](https://skillicons.dev/icons?i=py) Via pip 

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

Optional features require additional packages that are not needed for the base installation:

```bash
# ChEBI lookups
pip install "refineGEMs[chebi]"

# SBO label lookup via OLS
pip install "refineGEMs[ols]"

# SBO annotation
pip install "refineGEMs[sbo]"

# install all optional dependencies
pip install "refineGEMs[optional]"
```

> [!CAUTION]
> Some connected tools are optional and currently need to be installed directly from GitHub before using the
> corresponding refineGEMs workflow step. If they are missing, refineGEMs reports the missing dependency and skips
> the affected optional step where possible.
>
> ```bash
> # For MCC, until hot fix is merged into main:
> pip install "masschargecuration@git+https://github.com/Biomathsys/MassChargeCuration"
>
> # For BOFdat, our fork with hot fix(es):
> pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"
>
> # ModelPolisher client:
> pip install "model-polisher@git+https://github.com/draeger-lab/MPClient"
>
> ```

### ![docker](https://skillicons.dev/icons?i=docker) Via Docker

``refineGEMs`` can also be used via Docker. To build the docker image, firtsly clone the repository:

```bash
   git clone "https://github.com/draeger-lab/refinegems.git"
```

Then change into the directory and build the image:

```bash
   cd refinegems \
   docker build -t refinegems .
```
> [!NOTE]
> To provide the input files and retrieve the output files mount one folder as workspace folder to the Docker image with `-v`.

The default command executed by the image is ``refinegems -h`` and provides the help information for the CLI of 
``refineGEMs``.

```bash
   docker run refinegems -h
```

To use the image interactively and open a bash shell, run the following command:

```bash
   docker run -it --entrypoint bash refinegems
```

To use the image for specific commands, you can simply use every of the CLI commands as entrypoint. 
For example, to curate a (draft) model, run:

```bash
   docker run --name <container_name> -v <user_folder>:/rg_cont refinegems analyse stats ./path/to/model.xml
```

## License

The refineGEMs source code is distributed under the MIT license. Bundled
third-party data, database identifiers, adapted code, and connected external
tools remain under their own licenses or terms; see
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for details.


## How to cite
When using `refineGEMs`, please cite the latest publication:

Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas Dräger. 
Genome-scale metabolic models consistently predict in vitro characteristics of Corynebacterium
striatum. Front. Bioinform., oct 2023. [doi:10.3389/fbinf.2023.1214074](https://doi.org/10.3389/fbinf.2023.1214074).

## Repositories using refineGEMs
- [C_striatum_GEMs](https://github.com/draeger-lab/C_striatum_GEMs)
- draeger-lab/Cacnes - `private`
- draeger-lab/Cgranulosum - `private`
- draeger-lab/Koxytoca - `private`
- draeger-lab/Mfortuitum - `private`
- draeger-lab/Scohnii - `private`
- draeger-lab/Shaemolyticus - `private`
- draeger-lab/Ssanguinis - `private`
