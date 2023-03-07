Filling gaps with refineGEMs
============================

There are two possibilities to use refineGEMs to fill gaps.

Manual gap filling
------------------

.. autofunction:: refinegems.curate.add_reactions_from_table
    :noindex:

Automated gap filling
---------------------

The `gapfill` module was created to enable an automatic way of filling gaps in a model via genes.
This module can be used either with KEGG were you only need the KEGG organism ID or with BioCyc or with both (Options: 'KEGG', 'BioCyc', 'KEGG+BioCyc'). 

    .. warning:: 
        Current restrictions:
            - Only bacteria supported as missing reactions are filtered by compartments ('c', 'e', 'p').
            - Only models where the organisms have entries in KEGG and/or BioCyc can be gap filled.

    1. For the gap filling using BioCyc you will need to create an  `account <https://biocyc.org/>`. 
    2. Then you will search for the strain of your organism. 
    3. Then you click on `Tools` and select `Special SmartTables`.
       There you can download the tables "All genes of <organism>" and "All reactions of <organism>".   
    4. For the gene to reaction mapping table:
    
            i. Remove all columns except 'Gene Name',
            ii. then `choose a transform` ('Reactions of gene'), 
            iii. then add `property` ('Accession-2') 
            iv. and after that delete the 'Gene Name' column.
            v. Afterwards choose 'Accession-2',
            vi. and use the filter function to delete all empty rows.
            vii. Then `Export to Spreadsheet File` and choose `frame ids`.
            
    5. For the reactions table: 
    
        i. Remove all columns except 'Reaction',
        ii. then `choose a transform`: 
        
            a. 'Reactants of reaction'
            b. 'Products of reaction'
            
        iii. and then `choose property`: 
        
            a. 'EC-Number'
            b. 'Reaction-Direction'
            c. 'Spontaneous?' (property)
            
        iv. Afterwards `Export to Spreadsheet File` and choose `frame ids`.
        
    6. For the metabolites table: 
    
        i. Use MetaCyc database -> all compounds
        ii. Remove all columns except 'Compound',
        iii. then choose `property`:
        
            a. 'Object ID'
            b. 'Chemical Formula'
            c. 'InChI-Key'
            d. database links -> 'ChEBI'
            
        iv. Afterwards `Export to Spreadsheet File` and choose `common names`.
    
    Run times: 
    
        'KEGG': ~ 2h
        'BioCyc': ~ 45mins - 1h
        'KEGG+BioCyc': ~ 3 - 4h
        
You can either use it with refineGEMs, just edit the config ...

or you can use it in your local script with

.. code:: python
    from refinegems import gapfill as rgg

    rgg.gapfill_analysis(<your model>, gapfill_params, new filename)
