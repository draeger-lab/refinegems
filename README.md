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
```
python main.py --help
```
will give you access to all possible flags that can be used with the script.

```
Usage: main.py [OPTIONS]

  main function to run the program

Options:
  -i, --input TEXT         Path to GEM to be investigated
  -o, --output TEXT        Determine if output file should be created, default
                           FALSE

  -k, --keggpathways TEXT  new filename/path for model with KEGG pathways
  --help                   Show this message and exit.
```
