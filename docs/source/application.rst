Usage as standalone application
================================

The script ``main.py`` can be used directly in the command line after
entering the virtual environment with ``pipenv shell``.

The ``config.yaml`` file contains defaults for all variables that need
to be set by the user.

.. code:: yaml

  Description: > 
    This file can be adapted to choose what refinegems should do.
    Note: For windows use \ instead of / for the paths


  General Setting: >
    Path to GEM to be investigated

  model: 'data/e_coli_core.xml' #
  entrez_email: '' # necessary to access NCBI API
  visualize: TRUE


  Settings for scripts that manipulate the model: >
    They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.

  ### Addition of KEGG Pathways as Groups ###
  keggpathways: FALSE
  kegg_path:  '' # path where to save model with KEGG Groups

  ### SBO-Term Annotation ###
  sboterms: FALSE
  sbo_path: '' # path where to save model with sbo terms

  ### Model polishing ### The database of the model identifiers needs to be specified with 'id_db'
  polish: FALSE
  id_db: 'BIGG' # Required! 
  # Possible identifiers, currently: BiGG & VMH
  # For other IDs the `polish` function in `polish.py` might need adjustment
  lab_strain: FALSE # Needs to be set to ensure that protein IDs get the 'bqbiol:isHomologTo' qualifier
                    # & to set the locus_tag to the ones obtained by the annotation
  protein_fasta: '' # Path to used CarveMe input file, if exists; Needs to be set for lab_strain: True
  polish_path: '' # path where to save the polished model

  ### Charge correction ###
  charge_corr: FALSE
  charge_path: ''
  charge_report_path: ''

  ### Manual Curation ###
  man_cur: FALSE
  man_cur_type: 'gapfill' # either 'gapfill' or 'metabs'
  man_cur_table: 'data/manual_curation.xlsx'
  man_cur_path: '' # path where to save modified model

  Settings for scripts that investigate the model: >
    These are only necessary if none of the scripts to manipulate the model are used.

  # Set the out path for the csv created by multiple and a possible xlsx for model investigation
  out_path: ''

  # Set the basis medium to simulate growth from
  growth_basis: 'minimal_uptake' # 'default_uptake' or 'minimal_uptake'

  # Settings if you want to compare multiple models
  multiple: FALSE
  multiple_paths: ['data/e_coli_core.xml', '', '']

  # media to simulate growth from, available: SNM, LB, M9, SMM, CGXII, RPMI, no duplicates allowed
  media: ['SNM3', 'RPMI', 'CGXlab', 'LB', 'M9', 'CGXII', 'CasA']

  # determine whether the memote score should be calculated, default: FALSE
  memote: FALSE

  # Determine if output file for single model should be created, default: cl
  # Filename is set as the models name
  output: xlsx #cl, xlsx, csv 

  ### ModelSEED comparison ###
  modelseed: FALSE # set to False if not needed

The repository structure has the following intention: 

* ``refinegems/`` contains all the functions needed in ``main.py`` 
* ``data/`` contains all tables that are used by different parts of the script as well as a toy model ``e_coli_core.xml`` 
* Instead of using the files given in ``data/``, you can use your own files and just change the paths in ``config.yaml``. Please be aware that some functions rely on input in a certain format so make sure to check the files given in the ``data/`` folder and use the same formatting. 
* ``databases/`` contains the ``sql`` file as well as the ``db`` file necessary for the SBOAnn script by Elisabeth Fritze as well as the modules `gapfill`, `growth` and `modelseed`.
* The ``setup.py`` and ``pyproject.toml`` enable creating a PyPi package called ``refinegems``.