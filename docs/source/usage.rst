Usage
======

Usage as standalone application
-------------------------------

.. warning:: 
    | ``main.py`` will be deprecated from version 2.0.0 onwards.*
    | For the main pipeline in main.py see the CarveMe-ModelPolisher-based (CMPB) workflow in `SPECIMEN <https://github.com/draeger-lab/SPECIMEN>`__.
    | Beware that SPECIMEN only runs with refineGEMs 2.0.0 (starting with developmental versions).


The script ``main.py`` can be used directly in the command line after
entering the virtual environment with ``pipenv shell`` or ``conda activate <EnvName>``.

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
    This file can be adapted to choose what refineGEMs should do.
    Note: For windows use \ instead of / for the paths
  
  
  General Setting: >
    Path to GEM to be investigated
  
  model: 'data/e_coli_core.xml' 
  # Set the out path for all analysis files
  out_path: ''
  
  Settings for scripts that investigate the model: >
    These are only necessary if none of the scripts to manipulate the model are used.
  
  # Set to TRUE if you want pngs that aid in model investigation, will be saved to a folder called 'visualization'
  visualize: TRUE
  
  # Set the basis medium to simulate growth from
  growth_basis: 'minimal_uptake' # 'default_uptake' or 'minimal_uptake'
  
  # Set to TRUE if you want to simulate anaerobic growth
  anaerobic_growth: FALSE
  
  # Settings if you want to compare multiple models
  multiple: FALSE
  multiple_paths: # enter as many paths as you need below
    - 'data/e_coli_core.xml'
    - ''
    - ''
  single: TRUE # set to False if you only want to work with the multiple models
  
  # media to simulate growth from, just comment the media you do not want with a #
  media: 
    - 'LB'
    - 'RPMI'
    - 'M9'
    - 'SNM3'
    - 'CGXII'
    - 'CasA'
    - 'Blood'
    - 'dGMM'
    - 'MP-AU'
  
  # Determine whether the biomass function should be checked & normalised
  biomass: TRUE
  
  # determine whether the memote score should be calculated, default: FALSE
  memote: FALSE
  
  # compare metabolites to the ModelSEED database
  modelseed: FALSE # set to False if not needed
  
  
  Settings for scripts that manipulate the model: >
    They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.
  
  model_out: '' # path and filename to where to save the modified model
  entrez_email: '' # necessary to access NCBI API
  
  ### Addition of KEGG Pathways as Groups ###
  keggpathways: FALSE
  
  ### SBO-Term Annotation ###
  sboterms: FALSE
  
  ### Model polishing ### The database of the model identifiers needs to be specified with 'id_db'
  polish: FALSE
  id_db: 'BIGG' # Required!
  # Possible identifiers, currently: BiGG & VMH
  # For other IDs the `polish` function in `polish.py` might need adjustment
  lab_strain: FALSE # Needs to be set to ensure that protein IDs get the 'bqbiol:isHomologTo' qualifier
                    # & to set the locus_tag to the ones obtained by the annotation
  protein_fasta: '' # Path to used CarveMe input file, if exists; Needs to be set for lab_strain: True
  
  ### Charge correction ###
  charge_corr: FALSE
  
  ### Manual Curation ###
  man_cur: FALSE
  man_cur_type: 'gapfill' # either 'gapfill' or 'metabs'
  man_cur_table: 'data/manual_curation.xlsx'
  
  ### Automatic gap filling ###
  # All parameters are required for all db_to_compare choices except:
  # - organismid which is only required for db_to_compare: 'KEGG'/'KEGG+BioCyc'
  # - and biocyc_files which is not required for 'KEGG'
  gap_analysis: FALSE
  gap_analysis_params:
    db_to_compare: 'KEGG'  # One of the choices KEGG|BioCyc|KEGG+BioCyc
    organismid: 'T05059'  # Needs to be specified for KEGG
    gff_file: 'data/cstr.gff'  # Path to RefSeq GFF file 
    biocyc_files: 
      - 'Path0'  # Path to TXT file containing a SmartTable from BioCyc with the columns 'Accession-2' 'Reaction of gene' (-)
      - 'Path1'  # Path to TXT file containing a SmartTable with all reaction relevant information (*)
      - 'Path2'  # Path to TXT file containing a SmartTable with all metabolite relevant information (+)
      - 'Path3'  # Path to protein FASTA file used as input for CarveMe (Needed to get the protein IDs from the locus tags)
  # (-) If the organism is not in BioCyc retrieve a table mapping all reactions in BioCyc to the corresponding sequence
  # (*) 'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
  # (+) 'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
  gapfill_model: FALSE
  gap_analysis_file: 'Path to Excel file with which gaps in model should be filled'
  # Either obtained by running gapfill_analysis/Created by hand with the same structure as the result file from gapfill_analysis
  # Example Excel file to fill in by hand: data/modelName_gapfill_analysis_date_example.xlsx

The repository structure has the following intention: 

* ``refinegems/`` contains all the functions needed in ``main.py`` 
* ``data/`` contains all example tables that can be used as input for the curation scripts as well as the ``media_db.csv`` and a toy model ``e_coli_core.xml`` 
* Instead of using the files given in ``data/``, you can use your own files and just change the paths in ``config.yaml``. Please be aware that some functions rely on input in a certain format so make sure to check the files given in the ``data/`` folder and use the same formatting. 
* ``refinegems/databases/`` contains the SQL Schema file for the media and ``sboann``-related tables as well as the ready-to-use database file necessary for the SBOAnn script by Elisabeth Fritze as well as the modules ``gapfill``, ``growth`` and ``modelseed``.
* The ``setup.py`` and ``pyproject.toml`` enable creating a PyPI package called ``refineGEMs``.


Usage as python module
----------------------

.. warning:: 
    | *Function calls will be deprecated from version 2.0.0 onwards.*
    | Due to massive restructuring and extension, all functions have been moved into new modules.
    | If you used any functions from refineGEMs until now and want to use version 2.0.0 in the future, make sure to check your code after the version change.

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
.. autofunction:: refinegems.biomass.check_normalise_biomass
    :noindex:

The modules ``io``, ``cvterms`` and ``investigate`` provide functions that can be used by themselves.

``io`` 
^^^^^^
Provides a couple of helper functions to load models, databases and parse gfffiles

.. automodule:: refinegems.io
    :members:
    :noindex:

``investigate``
^^^^^^^^^^^^^^^
Provides a couple of functions to get parameters of your model

.. automodule:: refinegems.investigate
    :members:
    :noindex:

``cvterms``
^^^^^^^^^^^
Provides functions to work with cvterms

.. automodule:: refinegems.cvterms
    :members:
    :noindex:
