Welcome to refinegems!
======================================
``refinegems`` is a python package inteded
to help with the curation of genome-scale metabolic models (GEMS).

Overview
--------

Currently ``refinegems`` can be used for the investigation of a GEM, it can complete the following tasks:

* loading GEMS with ``cobrapy`` and ``libSBML``
* report number of metabolites, reactions and genes
* report orphaned, deadends and disconnected metabolites
* report mass and charge unbalanced reactions
* report `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ score
* compare the genes present in the model to the genes found in the `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Database (Note: this requires a gff file of your organism and the KEGG identifier of your organism)
* compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the `ModelSEED <https://modelseed.org/>`__ Database

Other applications of ``refinegems`` include curation of a given model these include: 

* correction of a model created with `CarveMe <https://github.com/cdanielmachado/carveme>`__ v.1.5.1 (for example moving all relevant information from the notes to the annotation field) this includes automated annotation of NCBI genes to the GeneProtein section of the model
* addition of `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Pathways as Groups (using the `libSBML <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html>`__ Groups Plugin)
* SBO-Term annotation based on a script by Elisabeth Fritze
* updating the annotation of metabolites and extending the model with reactions (for the purpose of gapfilling) based on a table filled by the user ``data/manual_annotations.xlsx``, note that this only works when the structure of the given table is used


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   application
   usage
   modules

* :ref:`genindex`
* :ref:`search`
