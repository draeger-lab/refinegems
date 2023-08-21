Welcome to refineGEMs!
======================================
``refineGEMs`` is a Python package intended
to help with the curation of genome-scale metabolic models (GEMS).

.. hint:: For bug reports please write issues on the `GitHub page <https://github.com/draeger-lab/refinegems/issues>`__ or open a discussion `here <https://github.com/draeger-lab/refinegems/discussions>`__.

Overview
--------

Currently ``refineGEMs`` can be used for the investigation of a GEM, it can complete the following tasks:

* loading GEMS with ``COBRApy`` and ``libSBML``
* report number of metabolites, reactions and genes
* report orphaned, deadends and disconnected metabolites
* report mass and charge unbalanced reactions
* report `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ score
* compare the genes present in the model to the genes found in:
  * the `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Database (Note: This requires a GFF file of your organism and the KEGG identifier of your organism.)
  * Or the `BioCyc <https://biocyc.org/>`__ Database (Note: This requires that a database entry for your organism exists in BioCyc.)
* compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the `ModelSEED <https://modelseed.org/>`__ Database.

Other applications of ``refineGEMs`` to curate a given model include: 

* The correction of a model created with `CarveMe <https://github.com/cdanielmachado/carveme>`__ v.1.5.1 (for example moving all relevant information from the notes to the annotation field) this includes automated annotation of NCBI genes to the GeneProduct section of the model,
* The addition of `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Pathways as Groups (using the `libSBML <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html>`__ Groups Plugin),
* Updating the SBO-Term annotations based on SBOannotator\ :footcite:p:`Leonidou2023_sboann`,
* Updating the annotation of metabolites and extending the model with reactions (for the purpose of filling gaps) based on a table filled by the user ``data/manual_annotations.xlsx`` (Note: This only works when the structure of the given table is used.),
* And extending the model with all information surrounding reactions including the corresponding GeneProducts and metabolites by filling in the table ``data/modelName_gapfill_analysis_date_example.xlsx`` (Note: This also only works when the structure of the given Excel file is used).


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   main-modules
   in_silico_media_generation
   API access <api>
   pipeline
   Notes for developers <development>

* :ref:`genindex`
* :ref:`search`

.. footbibliography::
