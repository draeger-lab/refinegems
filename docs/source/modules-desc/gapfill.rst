Filling gaps with refineGEMs
============================

There are two possibilities to use refineGEMs to fill gaps.

Manual gap filling
------------------

See :py:func:`~refinegems.curation.curate.add_reactions_from_table`



Automated gap filling
---------------------

The :py:mod:`~refinegems.curation.gapfill` module was created to enable an automatic way of filling gaps in a model via genes.

.. warning:: 
    Current restrictions:
        - Only bacteria supported as missing reactions are filtered by compartments ('c', 'e', 'p').
        - Only models where the organisms have entries in KEGG and/or BioCyc can be gap filled.
    
Run times:
    * 'KEGG': ~ 1h 10mins
    * 'BioCyc': ~ 10mins
    * 'KEGG+BioCyc': ~ 1h 20mins

The module :py:mod:`~refinegems.curation.gapfill` can be used to:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
a. perform gap analysis: :py:func:`~refinegems.curation.gapfill.gap_analysis`
b. add genes, metabolites and reactions from an Excel table to a model: :py:func:`~refinegems.curation.gapfill.gapfill_model`
c. perform gap analysis and add the result directly to a model: :py:func:`refinegems.curation.gapfill.gapfill`

.. warning:: 
    To use the gap analysis and directly add the result to a model, currently, one of the options 'BioCyc' or 'KEGG+BioCyc' has to be selected.
    For all other options the usage of gap_analysis combined with gapfill_model will result in an error.

Relevant parameters
^^^^^^^^^^^^^^^^^^^
To perform the gap analysis the following parameters are relevant for the config.yaml file:
(See :ref:`Data acquisition from BioCyc` on how to obtain the files 1 to 3)

.. code:: yaml

    gap_analysis: TRUE
        gap_analysis_params:
          db_to_compare: 'One of the choices KEGG|BioCyc|KEGG+BioCyc'
          organismid: 'KEGG Organism ID' # Needs to be specified for KEGG
          gff_file: 'Path to RefSeq GFF file' # Needs to be specified for KEGG 
          biocyc_files:
            - 'File 1: Path to gene to reaction mapping table'
            - 'File 2: Path to reaction table'
            - 'File 3: Path to compounds table'
            - 'File 4: Path to protein FASTA file used as input for CarveMe'

To add genes, metabolites and reactions from an Excel table to a model the following parameters need to be set:
(The Excel file is either obtained by running gapfill_analysis or created by hand with the same structure as the result file from gapfill_analysis.
An example Excel file to fill in by hand can be found in the cloned repository under ``data/modelName_gapfill_analysis_date_example.xlsx``)

.. code:: yaml

    gapfill_model: TRUE
        gap_analysis_file: 'Path to Excel file with which gaps in model should be filled'

Data acquisition from BioCyc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. If you have no BioCyc account you will need to create one. See `BioCyc /> Create Free Account <https://biocyc.org/new-account.shtml>` to create an account. 
2. Then you need to search for the strain of your organism.
3. Within the database of your organism you need to click on `Tools` in the menu bar and select `Special SmartTables` under `SmartTables`.
   There you need to make an adjustable copy of each of the tables "All genes of <organism>" and "All reactions of <organism>".   
4. **For the gene to reaction mapping table:**

        i. Remove all columns except 'Gene Name' from the "All genes of <organism>" table,
        ii. then click `choose a transform` and select 'Reactions of gene', 
        iii. then add the `property` 'Accession-2'
        iv. and delete the 'Gene Name' column.
        v. After that select the column 'Accession-2' and use the filter function in the box on the right side of the page to delete all empty rows.
        vi. Finally, click `Export to Spreadsheet File` from the box on the right side and choose `Frame IDs`.
        
5. **For the reactions table:** 

    i. Remove all columns except 'Reaction' from the "All reactions of <organism>" table,
    ii. then click `choose a transform`: 
    
        a. select 'Reactants of reaction',
        b. then select 'Products of reaction'
        
    iii. and then choose the `property`: 
    
        a. 'EC-Number',
        b. then 'Reaction-Direction',
        c. and then 'Spontaneous?'.
        
    iv. Finally, click `Export to Spreadsheet File` in the box on the right side and choose `Frame IDs`.
    
6. **For the metabolites table:** 

    i. Use the MetaCyc database to get the table "All compounds of MetaCyc".
    ii. Remove all columns except 'Compound',
    iii. then choose the `property`:
    
        a. 'Object ID',
        b. then 'Chemical Formula',
        c. then 'InChI-Key',
        d. and then 'database links' > 'ChEBI'.
        
    iv. Finally, click `Export to Spreadsheet File` in the box on the right side and choose `common names`.