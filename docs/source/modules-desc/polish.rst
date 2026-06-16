Polishing a draft model
=======================

The polishing workflow is intended for draft models created by automatic reconstruction
pipelines. It collects several recurring clean-up steps in
:py:func:`~refinegems.curation.curate.polish_model`, for example extending annotations,
cleaning notes fields, fixing qualifier patterns and adding model-level defaults that
are commonly missing in draft reconstructions.

The workflow can be run from the command line with:

.. code-block:: bash

    refinegems curate model MODEL.xml --outdir polished_model

For models with BiGG or VMH identifiers, the workflow can use the model entity IDs as
the main identifier source:

.. code-block:: bash

    refinegems curate model MODEL.xml --id-db BIGG --outdir polished_model

Gene product mappings
---------------------

Several polishing steps require a mapping between model gene products and external
identifiers such as RefSeq, NCBI protein identifiers or locus tags. The mapping can be
provided directly:

.. code-block:: bash

    refinegems curate model MODEL.xml --mapping-tbl-file geneproduct_mapping.tsv

Alternatively, the mapping table can be generated first from the model and one or more
GFF files:

.. code-block:: bash

    refinegems setup geneproduct-mapping-table MODEL.xml --gff-paths annotation.gff --outdir mapping

The generated table can then be checked manually and passed to ``curate model``.

For example, a small mapping table can be generated from the bundled *E. coli*
core model and an annotation file in GFF format:

.. code-block:: bash

    refinegems setup geneproduct-mapping-table e_coli_core.xml \
        --gff-paths ecoli_k12_mg1655_two_cds.gff \
        --contains-locus-tags True \
        --outdir mapping

The GFF file used here was shortened to two CDS entries from the NCBI RefSeq
annotation for *Escherichia coli* K-12 substr. MG1655
(`GCF_000005845.2_ASM584v2 <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz>`__,
chromosome ``NC_000913.3``):

.. code-block:: text

    ##gff-version 3
    #!genome-build ASM584v2
    #!genome-build-accession NCBI_Assembly:GCF_000005845.2
    ##sequence-region NC_000913.3 1 4641652
    NC_000913.3  RefSeq  CDS  372921   373871   .  +  0  ID=cds-NP_414885.1;Dbxref=GenBank:NP_414885.1,GeneID:945008;gene=mhpF;locus_tag=b0351;product=acetaldehyde dehydrogenase (acetylating) MhpF;protein_id=NP_414885.1
    NC_000913.3  RefSeq  CDS  1295446  1298121  .  -  0  ID=cds-NP_415757.1;Dbxref=GenBank:NP_415757.1,GeneID:945837;gene=adhE;locus_tag=b1241;product=fused acetaldehyde-CoA dehydrogenase and Fe-dependent alcohol dehydrogenasealdehyde/alcohol dehydrogenase AdhE;protein_id=NP_415757.1

The resulting table contains the model gene product IDs, matching locus tags,
classified database identifiers and gene product names. Rows without a match in
the shortened GFF are omitted from this example:

.. code-block:: text

    model_id,locus_tag,NCBI,REFSEQ,name
    G_b0351,b0351,b0351,NP_414885.1,acetaldehyde dehydrogenase (acetylating) MhpF
    G_b1241,b1241,b1241,NP_415757.1,fused acetaldehyde-CoA dehydrogenase and Fe-dependent alcohol dehydrogenasealdehyde/alcohol dehydrogenase AdhE

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:

    1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
        If there is no available data for the modeled organism in any database these identifiers can be added with
        the ``PGAB`` pipeline described in ``SPECIMEN`` From genome sequence to draft model` before draft model creation.
    2. Input of a FASTA file containing header lines similar to:
        >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
        Of the description part in the header line only locus_tag, protein and protein_id are important for ``cv_ncbiprotein``/ ``polish``.
