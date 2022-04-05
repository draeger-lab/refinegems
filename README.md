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
- correction of a model created with CarveMe v.1.5.1 (for example moving all relevant information from the notes to the annotation field)
- addition of [KEGG](https://www.genome.jp/kegg/kegg1.html) Pathways as Groups (using the [libSBML](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html) Groups Plugin)
- SBO-Term annotation based on a script by Elisabeth Fritze

## Installation

It is recommended to install all required packages in a `pipenv`. 
```
pip install pipenv
```
The `pipenv` package can also be installed via Anaconda (recommended if you are a Windows user).

All dependencies can the be installed by running 
```
pipenv install
```
in the local directory of this repository. All packages will be installed to pipenv called something like the repository name.

The command
```
pipenv shell
```
will initiate a session in the environment.

You should be all set now.

## Further Dependencies

If you want to use the SBO terms annotation part you will need to install PostgreSQL to your machine.

Make sure that you will be able to use psql from the command line. 

The database containing the BiGG ID and EC number mappings
```
sbo/create_dbs.sql
```
must be imported to a local PostgreSQL database to a selected user. 

You can use the following command, run it from this directory:
```
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
The script `main.py` can be used directly in the command line. 

The `config.yaml` file contains defaults for all variables that need to be set by the user.

```
# Path to GEM to be investigated
model: 'models/CStr_20210518.xml' 

# new filename/path for model with KEGG pathways, default: null means KEGG pathways will not be added
keggpathways: null 

# Requirements: PostgreSQL
sboterms: FALSE
database_user: postgres
database_name: sbo_ann
new_filename: 'cstr_sboann.xml'

### The following inputs are only necessary if neither kegg nor sbo is active ###

# Path to database which contains metabolites present in different media
media_db: 'media/media_db.csv' 

# media to simulate growth from, available: SNM, LB, M9, SMM, CGXII, RPMI
media: ['SNM3', 'LB', 'M9', 'SMM', 'CGXII', 'RPMI']

# determine whether the memote score should be calculated, default: FALSE
memote: FALSE

# Determine if output file should be created, default: cl
# Filename is set as the models name
output: cl #xlsx, csv 

### Gene comparison ###
# set to False if not needed
genecomp: TRUE
# the following is only relevant when turned on
organismid: 'T05059' # C. striatum
gff_file: 'genecomp/cstr.gff' # C. striatum
biggreactions: 'genecomp/bigg_models_reactions.txt'
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