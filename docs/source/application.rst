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
   entrez_email: 'famke.baeuerle@student.uni-tuebingen.de'


   Settings for scripts that manipulate the model: >
     They are all split into the ON / OFF switch (TRUE / FALSE) and additional settings like a path to where the new model should be saved.

   ### Addition of KEGG Pathways as Groups ###
   keggpathways: FALSE
   kegg_path: '../Nextcloud/master_thesis/models/Cstr_17_kegg.xml' # path where to save model with KEGG Groups

   ### SBO-Term Annotation (requires PostgreSQL) ###
   sboterms: FALSE
   sbo_path: '../Nextcloud/master_thesis/models/Cstr_17_sbo.xml' # path where to save model with sbo terms

   ### CarveMe polishing ###
   polish: FALSE
   id_db: 'BIGG' # Required!
   # Possible identifiers, currently: BiGG & VMH
   # For other IDs the `polish` function in `polish.py` might need adjustment
   lab_strain: FALSE # Needs to be set to ensure that protein IDs get the 'bqbiol:isHomologTo' qualifier
                  # & to set the locus_tag to the ones obtained by the annotation
   protein_fasta: '' # Path to used CarveMe input file, if exists; Needs to be set for lab_strain: True
   polish_path: '../Nextcloud/master_thesis/models/Cstr_17_genes.xml' # path where to save the polished model

   ### Charge correction ###
   charge_corr: FALSE
   charge_path: '../Nextcloud/master_thesis/models/Cstr_16_charges.xml'


   Settings for scripts that investigate the model: >
     These are only necessary if none of the scripts to manipulate the model are used.

   # media to simulate growth from, available: SNM, LB, M9, SMM, CGXII, RPMI
   media: ['SNM3', 'LB', 'M9']

   # determine whether the memote score should be calculated, default: FALSE
   memote: FALSE

   # Determine if output file should be created, default: cl
   # Filename is set as the models name
   output: xlsx #cl, xlsx, csv 

   ### Gapfill ###
   # All parameters are required for all db_to_compare choices except:
   # - organismid which is only required for db_to_compare: 'KEGG'/'KEGG+BioCyc'
   # - and biocyc_tables which is not required for 'KEGG'
   gapfill_analysis: FALSE
   gapfill_analysis_params:
     db_to_compare: 'KEGG+BioCyc'  # One of the choices KEGG|BioCyc|GFF|KEGG+BioCyc
     organismid: 'Three-letter KEGG OrganismID'  # Needs to be specified for KEGG
     gff_file: 'Path to RefSeq GFF file'  # Obtainable via the NCBI Assembly Accession page for the organism 
     biocyc_files: 
       - 'Path0'  # Path to TXT file containing a SmartTable from BioCyc with the columns 'Accession-2' 'Reaction of gene' (-)
       - 'Path1'  # Path to TXT file containing a SmartTable with all reaction relevant information (*)
       - 'Path2'  # Path to TXT file containing a SmartTable with all metabolite relevant information (+)
       - 'Path3'  # Path to protein FASTA file used as input for CarveMe (Needed to get the protein IDs from the locus tags)
   # (-) If the organism is not in BioCyc retrieve a table mapping all reactions in BioCyc to the corresponding sequence
   # (*) 'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
   # (+) 'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
   # For all BioCyc files the order should be the same. If the organism does not occur in the BioCyc database the
   # complete tables for reactions can be used with the same columns.
   gapfill_model: FALSE
   gapfill_model_in: 'Path to Excel file with which gaps in model should be filled' 
   # Either obtained by running gapfill_analysis/Created by hand with the same structure as the result file from gapfill_analysis
   # Example Excel file to fill in by hand: data/modelName_gapfill_analysis_date_example.xlsx
   gapfill_model_out: 'Path where to save gap filled model'

   ### ModelSEED comparison ###
   modelseed: TRUE # set to False if not needed

The repository structure has the following intention: 

* ``refinegems/`` contains all the functions needed in ``main.py`` 
* ``data/`` contains all tables that are used by different parts of the script as well as a toy model ``e_coli_core.xml`` 
* Instead of using the files given in ``data/``, you can use your own files and just change the paths in ``config.yaml``. Please be aware that some functions rely on input in a certain format so make sure to check the files given in the ``data/`` folder and use the same formatting. 
* ``databases/`` contains the ``sql`` file as well as the ``db`` file necessary for the SBOAnn script by Elisabeth Fritze as well as the modules `gapfill`, `growth` and `modelseed`.
* The ``setup.py`` and ``pyproject.toml`` enable creating a PyPi package called ``refinegems``.