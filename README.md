[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# gem_curation_template
Template repository structure for a genome-scale metabolic model curation project.

## Overview
`refinegems` is a python package inteded to help with the curation of genome-scale metabolic models (GEMS).

Currently `refinegems` can be used for the investigation of a GEM, it can complete the following tasks:
- loading GEMS with `cobrapy` and `libSBML`
- report number of metabolites, reactions and genes
- report orphaned, deadends and disconnected metabolites
- report mass and charge unbalanced reactions
- report [Memote](https://memote.readthedocs.io/en/latest/index.html) score
- compare the genes present in the model to the genes found in the [KEGG](https://www.genome.jp/kegg/kegg1.html) Database (Note: this requires a gff file of your organism and the KEGG identifier of your organism)
- compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the [ModelSEED](https://modelseed.org/) Database

Other applications of `refinegems` include curation of a given model these include:
- correction of a model created with [CarveMe](https://github.com/cdanielmachado/carveme) v.1.5.1 (for example moving all relevant information from the notes to the annotation field)
- addition of [KEGG](https://www.genome.jp/kegg/kegg1.html) Pathways as Groups (using the [libSBML](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html) Groups Plugin)
- SBO-Term annotation based on a script by Elisabeth Fritze

## Installation

`refinegems` is distributed via this github repository, all dependencies are denoted in a [pipenv](https://pipenv.pypa.io/en/latest/). You will need to install `pipenv` first. To install `refinegems` locally complete the following steps:

```bash
# install pipenv using pip
pip install pipenv

# clone or pull the latest source code
git clone https://github.com/draeger-lab/gem_curation_template.git
cd gem_curation_template

# install all dependencies
pipenv install

# initiate a session in the virtual environment
pipenv shell

```
The `pipenv` package can also be installed via Anaconda (recommended if you are a Windows user).

## Further Dependencies

If you want to use the SBO terms annotation part you will need to install PostgreSQL to your machine.

Make sure that you will be able to use psql from the command line. 

The database containing the BiGG ID and EC number mappings
```
sbo/create_dbs.sql
```
must be imported to a local PostgreSQL database to a selected user. 

You can use the following command, run it from this directory:
```bash
psql -U {your postgres username} -h localhost -d {your database name} < sbo/create_dbs.sql 
```

If you are a Windows user you will need to use a different command:
Enter into the psql shell by typing `psql`, then create the database with
```
CREATE DATABASE sbo_ann;
```
Afterwards load the database with
```
psql.exe -U postgres -d sbo_ann -f sbo\create_dbs.sql
```

## Usage
The script `main.py` can be used directly in the command line after entering the virtual environment with `pipenv shell`.

The `config.yaml` file contains defaults for all variables that need to be set by the user. 

```yaml
Description: > 
  This file can be adapted to choose what refinegems should do.
  Note: For windows use \ instead of / for the paths

---
General Setting: >
  Path to GEM to be investigated

model: 'models/cstr_ma.xml' #'models/CStr_20210518.xml' 

---
Settings for scripts that manipulate the model: >
  They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.

### Addition of KEGG Pathways as Groups ###
keggpathways: FALSE
kegg_path: ''

### SBO-Term Annotation (requires PostgreSQL) ###
sboterms: FALSE
{database_user: postgres, database_name: sbo_ann}
sbo_path: '../Nextcloud/master_thesis/models/Cstr_17_sbo.xml' # path where to save model with sbo terms

### CarveMe polishing ###
polish_carveme: FALSE
polish_path: '../Nextcloud/master_thesis/models/Cstr_17_clean.xml' # path where to save polished model

### Charge correction ###
charge_corr: FALSE
charge_path: '../Nextcloud/master_thesis/models/Cstr_17_char.xml'

---
Settings for scripts that investigate the model: >
  These are only necessary if none of the scripts to manipulate the model are used.

# Path to database which contains metabolites present in different media
media_db: 'media/media_db.csv' 

# media to simulate growth from, available: SNM, LB, M9, SMM, CGXII, RPMI
media: ['SNM3', 'LB', 'M9', 'SMM', 'CGXII', 'RPMI']

# determine whether the memote score should be calculated, default: FALSE
memote: FALSE

# Determine if output file should be created, default: cl
# Filename is set as the models name
output: xlsx #cl, xlsx, csv 

### Gene comparison ###
genecomp: FALSE # set to False if not needed
# the following is only relevant when turned on
organismid: 'T05059' # C. striatum
gff_file: 'genecomp/cstr.gff' # C. striatum
biggreactions: 'genecomp/bigg_models_reactions.txt'

### ModelSEED comparison ###
modelseed: FALSE # set to False if not needed
modelseedpath: 'modelseed/modelseed_compounds.tsv'
```

## Troubleshooting

* If you get `ImportError: DLL load failed while importing _sqlite3` when running main.py. Locate the `sqlite3.dll` file on you machine and add it to PATH.

* If you use python 3.8 it everything should work, just edit the `Pipfile` entry to `python_version = "3.8"` before running `pipenv install`.

* If you can't use `psql`from the command line, a common issue is that its not added to PATH:
```
locate psql | grep /bin
export PATH={Output from the line above with /bin as line end}:$PATH
```
* If you are a Windows user you will want to locate the installation manually. The path should look like something like this
`C:\Program Files\PostgreSQL\13\lib`.

* If you run into `psycopg2.OperationalError: fe_sendauth: no password supplied`: Change `scram-sha256`to `trust` in your file `pg_hba.conf` (located probably in `C:\Program Files\PostgreSQL\13\data`) 