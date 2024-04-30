Welcome to refineGEMs!
======================================
``refineGEMs`` is a Python package containing a collection of tools for the curation of genome-scale metabolic models (GEMS).

Overview
--------

Currently the ``refineGEMs`` toolbox includes the following features:

* loading GEMs with ``COBRApy`` and ``libSBML``
* report the statistics of a model, including 
  * number of metabolites, reactions and genes
  * number of orphaned, deadends and disconnected metabolites
  * number of mass and charge unbalanced reactions
* compare the genes present in the model to the genes found in:
  * the `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Database (Note: This requires the GFF file and the KEGG identifier of your organism.)
  * Or the `BioCyc <https://biocyc.org/>`__ Database (Note: This requires that a database entry for your organism exists in BioCyc.)
* compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the `ModelSEED <https://modelseed.org/>`__ Database.
* finding and resolving duplicates
* identifying and resolving EGCs (energy generating cycles)
* simulate growth under different conditions with an in.build media database and greatly customisable configurations
* and much more ...

 ``refineGEMs`` also provides access points for other tools and manual curation: 

* report `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ score or report
* Updating the SBO-Term annotations based on SBOannotator\ :footcite:p:`Leonidou2023_sboann`,
* The correction of a model created with `CarveMe <https://github.com/cdanielmachado/carveme>`__ v1.5.1 or v1.5.2 (for example moving all relevant information from the notes to the annotation field or automatically annotating the GeneProduct section of the model with the respective NCBI gene/protein identifiers from the GeneProduct identifiers),
* The addition of `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Pathways as Groups (using the `libSBML <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html>`__ Groups Plugin),
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
   database/intro
   Content of refineGEMs <api>
   pipeline
   Help and FAQ <help>
   Notes for developers <development>

* :ref:`genindex`
* :ref:`search`

.. footbibliography::
