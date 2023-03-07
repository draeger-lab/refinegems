Usage
======

Usage as standalone application
-------------------------------

The script ``main.py`` can be used directly in the command line after
entering the virtual environment with ``pipenv shell``.

The ``config.yaml`` file contains defaults for all variables that need
to be set by the user.

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:
        
      1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
          If there is no available data for the modeled organism in any database these identifiers can be added with 
          the pipeline described in :ref:`Pipeline: From genome sequence to draft model` before draft model creation.  
      2. Input of a FASTA file containing header lines similar to:
          >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
          Of the description part in the header line only locus_tag, protein and protein_id are important for ``polish``.
    

.. code:: yaml

  Description: > 
    This file can be adapted to choose what refinegems should do.
    Note: For windows use \ instead of / for the paths


  General Setting: >
    Path to GEM to be investigated

  model: 'data/e_coli_core.xml' 

  Settings for scripts that investigate the model: >
    These are only necessary if none of the scripts to manipulate the model are used.

  # Set to TRUE if you want pngs that aid in model inevstigation, will be saved to a folder called 'visualization'
  visualize: TRUE

  # Set the out path for the csv created by multiple and a possible xlsx for model investigation
  out_path: ''

  # Set the basis medium to simulate growth from
  growth_basis: 'minimal_uptake' # 'default_uptake' or 'minimal_uptake'

  # Settings if you want to compare multiple models
  multiple: FALSE
  multiple_paths: # enter as many paths as you need below
    - 'data/e_coli_core.xml'
    - ''
    - ''

  # media to simulate growth from, just comment the media you do not want with a #
  media: 
    - 'SNM3'
    - 'RPMI'
    - 'CGXlab'
    - 'LB'
    - 'M9'
    - 'CGXII'
    - 'CasA'

  # determine whether the memote score should be calculated, default: FALSE
  memote: FALSE

  # Determine if output file for single model should be created, default: cl
  # Filename is set as the models name
  output: xlsx #cl, xlsx, csv 

  # compare metabolites to the ModelSEED database
  modelseed: FALSE # set to False if not needed


  Settings for scripts that manipulate the model: >
    They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.

  entrez_email: '' # necessary to access NCBI API

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

The repository structure has the following intention: 

* ``refinegems/`` contains all the functions needed in ``main.py`` 
* ``data/`` contains all tables that are used by different parts of the script as well as a toy model ``e_coli_core.xml`` 
* Instead of using the files given in ``data/``, you can use your own files and just change the paths in ``config.yaml``. Please be aware that some functions rely on input in a certain format so make sure to check the files given in the ``data/`` folder and use the same formatting. 
* ``databases/`` contains the ``sql`` file as well as the ``db`` file necessary for the SBOAnn script by Elisabeth Fritze as well as the modules ``gapfill``, ``growth`` and ``modelseed``.
* The ``setup.py`` and ``pyproject.toml`` enable creating a PyPi package called ``refinegems``.


Usage as python module
----------------------

See :doc:`examples <modules/examples>` to learn how to use refineGEMs.
Note that at this time most of the modules only make sense when you use the respective main functions:

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:
        
    1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
        If there is no available data for the modeled organism in any database these identifiers can be added with 
        the pipeline described in :ref:`Pipeline: From genome sequence to draft model` before draft model creation.
        
    2. Input of a FASTA file containing header lines similar to:
        >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
        Of the description part in the header line only locus_tag, protein and protein_id are important for ``polish``.


.. autofunction:: refinegems.polish.polish
    :noindex:
.. autofunction:: refinegems.modelseed.compare_to_modelseed
    :noindex:
.. autofunction:: refinegems.pathways.kegg_pathways
    :noindex:
.. autofunction:: refinegems.sboann.sbo_annotation
    :noindex:

The modules ``io``, ``cvterms`` and ``investigate`` provide functions that can be used by themselves.

``io`` 
^^^^^^
provides a couple of helper functions to load models, databases and parse gfffiles

.. automodule:: refinegems.io
    :members:
    :noindex:

``investigate``
^^^^^^^^^^^^^^^
provides a couple of functions to get parameters of your model

.. automodule:: refinegems.investigate
    :members:
    :noindex:

``cvterms``
^^^^^^^^^^^
provides functions to work with cvterms

.. automodule:: refinegems.cvterms
    :members:
    :noindex:
