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

### Usage
The script `main.py` can be used directly in the command line. 

The `config.yaml` file contains defaults for all variables that need to be set by the user. The user can change the path to the model, the medium etc.
