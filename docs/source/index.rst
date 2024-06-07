Welcome to refineGEMs!
======================================
``refineGEMs`` is a Python-based toolbox for the curation and analysis of genome-scale metabolic models (GEMS).

Overview
--------

Currently the ``refineGEMs`` toolbox includes the following features:

* Loading GEMs with ``COBRApy`` and ``libSBML``
* Report the statistics of a model, including 

  * Number of metabolites, reactions and genes
  * Number of orphaned, deadends and disconnected metabolites
  * Number of mass and charge unbalanced reactions

* Compare the genes present in the model to the genes found in:

  * The `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Database (Note: This requires the RefSeq GFF file and the KEGG identifier of your organism.)
  * Or the `BioCyc <https://biocyc.org/>`__ Database (Note: This requires that a database entry for your organism exists in BioCyc.)

* Compare the charges and masses of the metabolites present in the model to the charges and masses denoted in the `ModelSEED <https://modelseed.org/>`__ Database.
* Finding and resolving duplicates
* Identifying and resolving EGCs (energy generating cycles)
* Simulate growth under different conditions with an in-build media database and greatly customisable configurations
* And much more ...

``refineGEMs`` also provides access points for other tools and manual curation: 

* Report `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ score or generate a complete `Memote <https://memote.readthedocs.io/en/latest/index.html>`__ report
* Updating the SBO-Term annotations based on `SBOannotator <https://github.com/draeger-lab/SBOannotator>`__\ :footcite:p:`Leonidou2023_sboann`,
* Balancing the masses and charges with `MassChargeCuration (MCC) <https://github.com/Biomathsys/MassChargeCuration>`__\ :footcite:p:`Finnem2023_mcc`,
* Correcting a biomass objective function with `BOFdat <https://github.com/jclachance/BOFdat>`__,
* The correction of a model created with `CarveMe <https://github.com/cdanielmachado/carveme>`__ v1.5.1, v1.5.2 or higher (for example moving all relevant information from the notes to the annotation field or automatically annotating the GeneProduct section of the model with the respective NCBI gene/protein identifiers from the GeneProduct identifiers),
* The addition of `KEGG <https://www.genome.jp/kegg/kegg1.html>`__ Pathways as Groups (using the `libSBML <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_groups_model_plugin.html>`__ Groups Plugin),
* Updating the annotation of metabolites and extending the model with reactions (for the purpose of filling gaps) based on a table filled by the user (Note: This only works when the structure of the `example manual curation Excel file <https://github.com/draeger-lab/refinegems/blob/dev-2/src/refinegems/example/example_inputs/manual_curation.xlsx>`__ is used.),
* And extending the model with all information surrounding reactions including the corresponding GeneProducts and metabolites by filling in the table (Note: This only works when the structure of the `example gapfill analysis Excel file <https://github.com/draeger-lab/refinegems/blob/dev-2/src/refinegems/example/example_inputs/modelName_gapfill_analysis_date_example.xlsx>`__ is used).

How to cite
-----------

When using ``refineGEMs``, please cite the latest publication:

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
