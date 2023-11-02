Welcome to refineGEMs!
======================================
``refineGEMs`` is a Python package intended
to help with the curation of genome-scale metabolic models (GEMS).

.. hint:: For bug reports please write issues on the `GitHub page <https://github.com/draeger-lab/refinegems/issues>`__ or open a discussion `here <https://github.com/draeger-lab/refinegems/discussions>`__.

Overview
--------

Currently ``refineGEMs`` can be used for the investigation of a GEM, it can complete the following tasks:

* loading GEMs with ``COBRApy`` and ``libSBML``
* report number of metabolites, reactions and genes
* report orphaned, deadends and disconnected metabolites
* report mass and charge unbalanced reactions
* report `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ score
* compare the genes present in the model to the genes found in:
  * the `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Database (Note: This requires the GFF file and the KEGG identifier of your organism.)
  * Or the `BioCyc <https://biocyc.org/>`__ Database (Note: This requires that a database entry for your organism exists in BioCyc.)
* compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the `ModelSEED <https://modelseed.org/>`__ Database.

Other applications of ``refineGEMs`` to curate a given model include: 

* The correction of a model created with `CarveMe <https://github.com/cdanielmachado/carveme>`__ v1.5.1 or v1.5.2 (for example moving all relevant information from the notes to the annotation field or automatically annotating the GeneProduct section of the model with the respective NCBI gene/protein identifiers from the GeneProduct identifiers),
* The addition of `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Pathways as Groups (using the `libSBML <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html>`__ Groups Plugin),
* Updating the SBO-Term annotations based on SBOannotator\ :footcite:p:`Leonidou2023_sboann`,
* Updating the annotation of metabolites and extending the model with reactions (for the purpose of filling gaps) based on a table filled by the user ``data/manual_annotations.xlsx`` (Note: This only works when the structure of the `example Excel file <https://github.com/draeger-lab/refinegems/blob/5eac900d9848b5ae5faf0055db72a986e7ba64e8/data/manual_curation.xlsx>`__ is used.),
* And extending the model with all information surrounding reactions including the corresponding GeneProducts and metabolites by filling in the table ``data/modelName_gapfill_analysis_date_example.xlsx`` (Note: This also only works when the structure of the `example Excel file <https://github.com/draeger-lab/refinegems/blob/5eac900d9848b5ae5faf0055db72a986e7ba64e8/data/modelName_gapfill_analysis_date_example.xlsx>`__ is used).

How to cite
-----------

When using refineGEMs, please cite the latest publication:

  Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas Dräger. 
  Genome-scale metabolic models consistently predict in vitro characteristics of Corynebacterium
  striatum. Front. Bioinform., oct 2023. `doi:10.3389/fbinf.2023.1214074 <https://doi.org/10.3389/fbinf.2023.1214074>`__.



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
