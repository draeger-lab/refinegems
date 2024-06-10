Help and FAQ
============

This page provides furth help if you run into problems with refineGEMs.

.. warning::

    With the release of version 2.0, refineGEMs will no longer work with Python ....

Help
----




Known Issues and Bugs
---------------------

.. hint:: 
    For bug reports please write issues on the `GitHub page <https://github.com/draeger-lab/refinegems/issues>`__ 
    or open a discussion `here <https://github.com/draeger-lab/refinegems/discussions>`__.

Pydantic
^^^^^^^^

Pydantic warning `underscore_attrs_are_private has been removed` has not - yet - caused any issues
however, the core of the problems (= what causes the warning) has yet to be identifies. 

FAQ
---

**When do I use ``lab_strain=True``?**

Using ``lab_strain=True`` has the following two requirements:
    
    1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
        If there is no available data for the modeled organism in any database these identifiers can be added with 
        the ``PGAB`` pipeline described in ``SPECIMEN`` before draft model creation.  
    2. Input of a FASTA file containing header lines similar to:
        >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
        Of the description part in the header line only locus_tag, protein and protein_id are important for :py:mod:`~refinegems.curation.polish`.
