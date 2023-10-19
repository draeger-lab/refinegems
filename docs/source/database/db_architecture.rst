*In silico* media database
==========================

This section describes the architecture of the database, as well as information on what and how the data is stored.

Entity-Relationship-Diagram
---------------------------

.. image:: ../images/ER_database.pdf
    :alt: Entity-Relationship-diagram for the database with examplary tables.

The diagram above shows the Entity-Relationship-diagram \(ER-diagram\) for the database, as well as examples of the tables that the final database consists of.
The ER-diagram shows the three entities medium, substance and database \(blue\), there relations \(orange\) and attributes \(green˜). 
Each medium includes multiple substances from different biological origins that can be absorbed by the model with the stated fluxes.
Since multiple media can include the same substance, the relationship is marked as n:m. 
Each substance is described by a name and a formula in addition to an unique identifier. 
Furthermore, the substances can be found in different databases and connected to the respective identifiers. 
Since this 1:n relation does not include further attributes, the tables was combinded with the database one, forming a total of four tables in which the data is stored.

Medium
------

The medium table contains an entry for each medium in the database, giving it an unique ID and storing its name and description.
The medium name should be a short abbreviation that is generally used to refer to this medium, e.g. LB.
The medium description should be a short sentence explaining its abbreviation, e.g. for LB the description "lysogengy broth" explains the abbreviation used as the name.

Medium_includes_substance \(m2s\)
---------------------------------

This table links the medium ID to the substances contained in the medium. Additionally, the flux for the possible exchange reaction is stored. 
In case the original publication the medium did not provide information for the flux, the attribite is left empty \(NULL\).
Furthermore, since the chemical substance in the medium might differ from the compund used to create the medium, the column origin provides the information about how and in which form the substance was originally added to the medium.

Substance 
---------

This table consists of the substamce ID, its name and its chemical formula. 
The chemical formula is saved without charge. 
In case of substances like vitamin B12, no substance formula is saved as these describe a group of different metabolites, 
that are treated as a single substance in some reactions.

Substaces are added using the following naming conventions:

1. Names are mostly based on CHEBI 
2. For amino acids, hormones and the like, use trivial name, otherwise use CHEBI IUPAC name with most used alias in square brackets 
3. For acids, use their salt form, e.g., sulfuric acid --> sulfate
4. | Metallic ions are entered using the following scheme: 
   | Name\(charge in roman numerals\) cat/anion \[chemical formula\(charge in arabic numerals\)\]
   | e.g. Calcium\(II\) cation \[Ca\(2+\)\]
5. | Halogenic ions are saved in a similar matter:
   | e.g. Chloride \[Cl\(-\)\]

.. note::
    When coverting the substances from the previous database format to this, well-annotated names were taken from the
    previous table as they were - meaning that they not necessarily follow the rules defined above. 
    These names might be subject to future change. 
    Any newly added substance should follow the namescheme above.
    For more information about the conversion, see the jupyter notebook **Remapping_substance_table.ipynb**

.. hint::
    | Before adding a new substance to the table, please check carefully if the substance already has an entry and if so, reuse it.
    | Help regarding mapping your infomration e.g. the BiGG to a formula or the CHEBI IUPAC name, please refer to the functions XXX and YYY.


Substance2db
------------

.. warning:: 
    This section is still incomplete and object to change.

The final table of the database contains the mapping of the substance identifiers to different database identifiers. 
To ensure no problems occurs in the case multiple databases have the same identifiers and help with searching for a certain kind of ID, the database type is stored in an additional column.
Currently, IDs from the following databases can be found in this table:

- BiGG
- SEED
- KEGG
- 

.. note::
    Be aware that it is possible for some substances to be missing IDs for certain databases in the list above as they either do not exist or have yet to be entered into the database.
