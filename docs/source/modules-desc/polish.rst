Polishing a CarveMe model
=========================

CarveMe version 1.5.1 leads to some irritations in the model, the scripts in 
:py:mod:`~refinegems.curation.polish` enable for example the addition of BiGG IDs 
to the annotations as well as a correct formatting of the annotations.

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:
        
    1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
        If there is no available data for the modeled organism in any database these identifiers can be added with 
        the pipeline described in :ref:`Pipeline: From genome sequence to draft model` before draft model creation.
    2. Input of a FASTA file containing header lines similar to:
        >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
        Of the description part in the header line only locus_tag, protein and protein_id are important for ``cv_ncbiprotein``/ ``polish``.
        