## gem_curation_template
Template repository structure for a genome-scale metabolic model curation project.

### Installation

It is recommended to install all required packages in a `pipenv`. 
```
pip install pipenv
```

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

### Further Dependencies

If you want to use the SBO terms annotation part you will need to install PostgreSQL to your machine.

Make sure that you will be able to use psql from the command line. A common issue is that its not added to PATH:
```
locate psql | grep /bin
export PATH={Output from the line above with /bin as line end}:$PATH
```

The database containing the BiGG ID and EC number mappings
```
sbo/create_dbs.sql
```
must be imported to a local PostgreSQL database to a selected user. 

You can use the following command, run it from this directory:
```
psql -U {your postgres username} -h localhost -d {your database name} < sbo/create_dbs.sql 
```

### Usage
The script `main.py` can be used directly in the command line. 

The `config.yaml` file contains defaults for all variables that need to be set by the user.

```
# Path to GEM to be investigated
model: 'models/CStr_20210518.xml' 

# new filename/path for model with KEGG pathways, default: null means KEGG pathways will not be added
keggpathways: null 

# Requirements: PostgreSQL
# needs database name and user: [database_user, database_name, 'new_filename.xml']
sboterms: null

### The following inputs are only necessary if neither kegg nor sbo is active

# Path to database which contains metabolites present in different media
media_db: 'media/media_db.csv' 

# media to simulate growth from, available: SNM, LB, M9, SMM
media: ['SNM3', 'LB', 'M9', 'SMM']

# determine whether the memote score should be calculated, default: FALSE
memote: FALSE

# Determine if output file should be created, default: [command_line, 0]
# for output into an excel file indicate the filename in this format: ['path/to/file/name.xlsx', 1]
# for output of the growth simulation into a csv file: ['path/to/file/name.csv', 2]
output: [command_line, 0]
```
