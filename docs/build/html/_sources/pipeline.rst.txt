Pipeline: From genome sequence to draft model
=============================================

Generating a model for an organism where no information on genes and proteins is obtainable via any database 
causes the problem that the model will not contain valid database identifiers for any GeneProduct. To resolve this issue the 
workflow in Figure :numref:`workflow` can be used.

1. First annotate the genome with NCBI's Prokaryotic Genome Annotation Pipeline (PGAP) to obtain the same FASTA format as used in NCBI.
2. Then use diamond with the ``nr`` database from NCBI and the obtained annotated FASTA file as input. Restrict the search to your organism's taxon if known and use the flag for taxonomy checking.
3. Check if any protein in the annotation FASTA file still has no database identifier.

    | -> YES: Rerun diamond without the taxonomy check and without the restriction for the organism's taxon.
    |
    | -> NO: Continue with step 4.

4. Add the diamond result to the annotated FASTA file.
5. Run e.g. ``CarveME``to obtain a draft model.
6. Check if in the model any GeneProducts without NCBI Protein or RefSeq identifiers occur.

    | -> YES: 
    |     i. Use individual BLAST searches for the remaining GeneProducts.
    |     ii. Add the results to the annotated and refined FASTA file.
    |     iii. Create again a draft model with the same program with the newly refined FASTA file.
    | 
    | -> NO: The draft model is done.

.. _workflow:
.. figure:: images/genome2draft.svg
  :alt: Workflow from genome sequence to a draft model

  Workflow from genome sequence to a draft model
