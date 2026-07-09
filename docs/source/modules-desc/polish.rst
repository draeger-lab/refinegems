Polishing a draft model
=======================

The polishing function is intended for draft models created by automatic reconstruction
workflows. It collects several recurring clean-up steps in
:py:func:`~refinegems.curation.curate.polish_model`, for example extending annotations,
cleaning notes fields, fixing qualifier patterns and adding model-level defaults that
are commonly missing in draft reconstructions. In the following mostly examples for 
the command line are shown. All options available via the command line are also 
available in the Python function.

The function can be run from the command line with defaults via

.. code-block:: bash

    refinegems curate model MODEL.xml --outdir polished_model

or directly used in Python via

.. code-block:: python

    from refinegems.curation.curate import polish_model

    polish_model("MODEL.xml", outdir="polished_model")

The function uses the model entity IDs. To specify the namespace of the identifiers use ``--id-db``. 
The default is ``BiGG``.

.. note::

    Currently, the implementation is only tested for models with identifiers from the BiGG Models database or VMH.

.. code-block:: bash

    refinegems curate model MODEL.xml --id-db BIGG --outdir polished_model

To see the other options simply call the help message of ``curate model``:

.. code-block:: bash

    refinegems curate model --help

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:

    1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
        If there is no available data for the modelled organism in any database these identifiers can be added with
        the `PGAB pipeline described in SPECIMEN <https://specimen.readthedocs.io/en/latest/pipeline_idea.html>`__ 
        before draft model creation.
    2. Input of a FASTA file containing header lines similar to:
        >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]

        Of the description part in the header line only locus_tag, protein and protein_id are important for ``cv_ncbiprotein``/ ``polish``.

Enhancing the gene product annotations
--------------------------------------

Within the ``polish_model`` function, additional annotations for the gene products of the provided model can be added.
This is done either via a mapping table for RefSeq, NCBI Protein IDs and locus tags or by querying the KEGG database for 
KEGG and UniProt IDs. The mapping table can be provided directly with the ``--mapping-tbl-file`` option to 
``curate model``:

.. code-block:: bash

    refinegems curate model MODEL.xml --mapping-tbl-file geneproduct_mapping.tsv

Alternatively, the mapping table can be generated from the model, one or more
GFF files (optional) and by querying the NCBI database (optional). This can be 
accessed via

.. code-block:: bash

    refinegems setup geneproduct-mapping-table MODEL.xml --gff-paths annotation.gff -email user@example.com --outdir mapping

or

.. code-block:: python

    from refinegems.utility.entities import generate_geneproduct_mapping_table

    get_gpid_mapping("MODEL.xml", gff_paths=["annotation.gff"], email="user@example.com", outpath="mapping")

The generated table can then be checked manually and passed to ``curate model`` or ``polish_model``.

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

.. hint::

    The mapping file does not need to be generated and provided beforehand to ``curate model``.
    If the mapping file is not provided, but an enhacement of the gene product annotations is desired, the mapping table 
    can be generated automatically by ``curate model``, if the corresponding flags are set.
    See the flags under `Parameters required when --mapping-tbl-file is not provided:` in the help message of ``curate model``.
