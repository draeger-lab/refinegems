Usage as python module
======================

See :doc:`examples <examples>` to learn how to use refineGEMs.
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
------
provides a couple of helper functions to load models, databases and parse gfffiles

.. automodule:: refinegems.io
    :members:
    :noindex:

``investigate``
---------------
provides a couple of functions to get parameters of your model

.. automodule:: refinegems.investigate
    :members:
    :noindex:

``cvterms``
-----------
provides functions to work with cvterms

.. automodule:: refinegems.cvterms
    :members:
    :noindex:
