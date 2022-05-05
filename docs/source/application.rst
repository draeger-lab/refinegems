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

   model: 'data/e_coli_core.xml' #'../C_striatum_GEMs/models/Cstr_14.xml'


   Settings for scripts that manipulate the model: >
     They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.

   ### Addition of KEGG Pathways as Groups ###
   keggpathways: FALSE
   kegg_path: '../Nextcloud/master_thesis/models/Cstr_17_kegg.xml' # path where to save model with KEGG Groups

   ### SBO-Term Annotation (requires PostgreSQL) ###
   sboterms: FALSE
   database_user: postgres
   database_name: sbo_ann
   sbo_path: '../Nextcloud/master_thesis/models/Cstr_17_sbo.xml' # path where to save model with sbo terms

   ### CarveMe polishing ###
   polish_carveme: FALSE
   polish_path: '../Nextcloud/master_thesis/models/Cstr_17_genes.xml' # path where to save the polished model
   entrez_email: 'famke.baeuerle@student.uni-tuebingen.de'

   ### Charge correction ###
   charge_corr: FALSE
   charge_path: '../Nextcloud/master_thesis/models/Cstr_16_charges.xml'


   Settings for scripts that investigate the model: >
     These are only necessary if none of the scripts to manipulate the model are used.

   # Path to database which contains metabolites present in different media
   media_db: 'data/media_db.csv' 

   # media to simulate growth from, available: SNM, LB, M9, SMM, CGXII, RPMI
   media: ['SNM3', 'LB', 'M9']

   # determine whether the memote score should be calculated, default: FALSE
   memote: FALSE

   # Determine if output file should be created, default: cl
   # Filename is set as the models name
   output: xlsx #cl, xlsx, csv 

   ### Gene comparison ###
   genecomp: FALSE # set to False if not needed
   # the following is only relevant when turned on
   organismid: 'T05059' # C. striatum
   gff_file: 'data/cstr.gff' # C. striatum
   biggreactions: 'data/bigg_models_reactions.txt'

   ### ModelSEED comparison ###
   modelseed: TRUE # set to False if not needed
   modelseedpath: 'data/modelseed_compounds.tsv'

The repository structure has the following intention: 

* ``refinegems/``
contains all the functions needed in ``main.py`` 
* ``data/`` contains
all tables that are used by different parts of the script as well as a
toy model ``e_coli_core.xml`` 
* Instead of using the files given in
``data/``, you can use your own files and just change the paths in
``config.yaml``. Please be aware that some functions rely on input in a
certain format so make sure to check the files given in the ``data/``
folder and use the same formatting. 
* ``sbo/`` contains the ``sql``
files necessary for the SBOAnn script by Elisabeth Fritze 
*
``setup.py`` and ``pyproject.toml`` enable creating a PyPi package
called ``refinegems``