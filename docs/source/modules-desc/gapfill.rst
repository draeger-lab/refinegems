Filling gaps with refineGEMs
============================

Finding and filling gaps in a genome-scale metabolic model is a frequently discussed and 
time-consuming part of the modelling process.

The :py:mod:`~refinegems.classes.gapfill` module of ``refineGEMs`` provides different 
flavors of :ref:`(Semi-)automated gap filling algorithmns`.

.. warning:: 

    The gap filling has undergone major restructering. Gap filling of version 
    > *2.0.0* behaves fundamentally different than the implementations in older versions. 

Currently, ``refineGEMs`` includes three ways of (semi-)automated gap filling:

- | :ref:`Gap filling with KEGG`:
  | If the organism to be modelled has a *KEGG organism ID*, this ID can be used to extract the genes and related enzymes from the KEGG database that are missing in the model and attempt to add them to it.

- | :ref:`Gap filling with BioCyc`:
  | If the organism to be modelled has an *entry in BioCyc*, this information can be compared to the model to add missing genes, reactions and more.

- | :ref:`Gap filling with a GFF (and SwissProt)`:
  | This algorithmn takes the protein GFF file of the organism and blasts the missing genes (products) against the SwissProt database to find homolgs, that can then be added to the model.

----

(Semi-)automated gap filling algorithmns
----------------------------------------

The idea behind these algorithmns is to reduce the amount of manual curation as much as 
possible without losing information along the way. Which algorithmn to choose mainly 
depends on the available information. Running multiple algorithmns subsequently is also 
possible.

The algorithmns have the same basic architecture:

1. Based on the available information, identify missing genes in the model (locus tag and protein (NCBI) accession number).
2. Map the missing genes to more information, e.g. the EC number, to identfy missing reactions corresponding to the gene products.
3. Use the collected information to try and fill the gaps in the model. 

    a. Add genes to gene production rules (GPRs) that have been mapped to existing reactions.
    b. Try to add reactions and metabolites to the model, if they can be successfully (completely) constructed in an automated manner.
    c. Add the genes and GPRs of the reactions added in the previous step.

During the steps above, if the information content if insufficient or the constructions of a model entity fails, 
the corresponding part is collected and saved for the user to enable faster manual curation afterwards (if desired).

All algorthmns have similar table outputs for the finding part, the filling part is therefore
same for all (see :py:meth:`~refinegems.classes.gapfill.GapFiller.fill_model`). 
Before adding any reactions to the model, it is checked, if it already exists to
reduce the possibility of duplicates being added.

Additionally, the following parameters allow the user to set up restrictions on 
which reactions and metabolites that should added to the model:

- ``exclude_dna``: If set to True, reactions containing the keyword ``DNA`` in their name are not added.
- ``exclude_rna``: If set to True, reactions containing the keyword ``RNA`` in their name are not added.
- ``formula_check``: Allow only reactions, whoose metabolites' formulas are of a certain quality.
    
    a. ``"none"``: No restriction
    b. ``"existence"``: Formula exists (exludes empty string and None/NaN values)
    c. ``"wildcard"``: Formula exists and does not contain the wildcard symbol ``"*"``
    d. ``"strict"``: Extends the previous option to also exclude formulas with a rest, denoted as ``"R"``

Gap filling with KEGG
^^^^^^^^^^^^^^^^^^^^^

| **Requirement:** KEGG organism ID 
| **Class:** :py:class:`~refinegems.classes.gapfill.KEGGapFiller`
| **Runtime estimation:** *to be determined*

To find the missing genes, the genes in the model are compared to the ones that can be
extracted from KEGG with the given organism ID. The comparison is based on the KEGG 
locus tags (format :code:`<kegg-organism-id>:<locus-tag>`). The IDs for the missing
genes are then used to retrieve the corresponding KEGG entry to extract information 
about related enzymes and reactions (via EC number and KEGG reaction ID). If a KEGG 
reaction ID is found, it can be directly used as a missing reactions, if an EC number s found, 
if is queries to KEGG to retrieve the reaction information corresponding to this 
number. 


Gap filling with BioCyc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Requirement:** BioCyc entry for the organism, access to BioCyc smart tables
| **Class:** :py:class:`~refinegems.classes.gapfill.BioCycGapFiller`
| **Runtime estimation:** *to be determined*

@TODO this fine?
@TODO Add description on how to get the tables!!!

If an organism has an entry for its metabolism in BioCyc, one can download two smart tables 
containing the available information about the genes (at least the columns ``Accession-2`` and 
``Reactions of gene``) and the reactions (at least the columns ``Reaction | Object ID | EC-Number | Spontanous?``).

These two tables, together with the GFF file are the required input for this gap filling algorithmn.
The missing genes are identfied by comparing the gene table ``Accession-2`` column to the model.
Subsequently, the missing genes are mapped back to the reactions to identify missing reactions.
The reactions are further mapped to MetaNetX and BiGG to obtain more reaction equations and 
information, since especially the metabolites are easier to construct using the other databases. 


Gap filling with a GFF (and SwissProt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Requirement:** Protein GFF (RefSeq or GenBank format)
| **Class:** :py:class:`~refinegems.classes.gapfill.GeneGapFiller`
| **Runtime estimation:** *to be determined*

In contrast to the other gap filling options, this one can be applied, if the organism has no database entry. 
Therefore, this gap filling algorithmn also works with newly discovered strains.

The idea is to extract the coding sequences of the organism from the GFF and map the corresponding
locus tags to the ones found in the model to identfy missing genes. Subsequently, the sequences of the 
missing genes are blasted against the SwissProt database to identfy homologs. The homologs are then mapped to
EC numbers (if possible). If the GFF already contains EC number information, these are extracted beforehand
to reduce the number of sequences, that need to be blasted. Additionally, the (NCBI) protein IDs 
can be searched in NCBI to extract information from there. This behaviour can be useful, if 
the input is a RefSeq GFF. It can be enabled by passing an e-mail address to the parameter :code:`mail` and 
setting :code:`check_NCBI` to `True` when running :py:meth:`~refinegems.classes.gapfill.GeneGapFiller.find_missing_reacs`. 
Finally, the EC numbers are mapped to different databases to find the
reactions, that should be added to the model. 

.. note:: 

    Currently, this way of gap filling is restricted to the SwissProt database, 
    however, there are plans to extend it to allow other/multiple databases 
    for the mapping.


How to run a GapFiller
----------------------

Due to the gap filling algorithmns having the same architecture, the function calls
for running them are basically the same, save for some parameters (will be denoted as ``<params>`` 
in the folowing code snippets.)

Firstly, the class instance for the chosen gapfiller, denoted by the place holder 
``<CHOSEN_GAPFILLER>``, must be initialed.

.. code-block:: python 
    :class: copyable
    
    from refinegems.classes.gapfill import <CHOSEN_GAPFILLER> # e.g. GeneGapFiller

    gapfiller = <CHOSEN_GAPFILLER>(<params>) 

The next step is to identfy the missing genes. Depending on the algorithm, some
additional parameters need to be added.

.. code-block:: python 
    :class: copyable

    # model = model loaded with libsbml
    gapfiller.find_missing_genes(model, <params>)    

Then, the missing reactions are identfied in a similar matter. The biggest differences
is, that this part relies on the model loaded with COBRApy, while the gene-finding part 
relies on the model loaded with libSBML. 

.. code-block:: python
    :class: copyable

    # cobramodel = model loaded with cobrapy
    gapfiller.find_missing_reacs(cobramodel, <params>)

Finally, the model can be extended with the collected information - as much as is automatic possible.

.. code-block:: python 
    :class: copyable

    # any_model = model loaded with either libsbml or cobrapy
    filled_model = gapfiller.fill_model(any_model, <params>)

To access information between steps or afterwards, the following attributes can be of interest:

    - :code:`gapfiller.missing_genes`: Table of currently missing and not further categorised genes.
    - :code:`gapfiller.missing_reacs`: Table of currently missing and not further categorised reactions.
    - :code:`gapfiller._statistics`: Dictionary of statistical values, e.g. number of added genes.
    - :code:`gapfiller.manual_curation`: Dictionary of tables containing information, that cannot be added automatically due to different reasons. Reason is denoted in the key.

    Some GapFillers also provide additional, for the corresponding algorithm specific, attributes.