#!/usr/bin/env python
"""The `gapfill` module can be used either with KEGG were you only need the KEGG organism ID or with BioCyc or with both (Options: 'KEGG', 'BioCyc', 'KEGG+BioCyc'). 

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
        
"""
import ast
from libsbml import Model
from libsbml import *
import refinegems.analysis_kegg as rga_kegg
import refinegems.analysis_biocyc as rga_biocyc
from refinegems.cvterms import add_cv_term_metabolites, add_cv_term_reactions 
from refinegems.entities import create_gp, create_species, create_reaction
import pandas as pd
import numpy as np
from typing import Union
from colorama import init as colorama_init
from colorama import Fore

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"

'''Skeleton for functions that could be used for a lab strain/organism which is in no database contained
def get_genes_from_gff():
    pass


def get_related_metabs_reactions_blast():
    pass


def gff_gene_comp():
    pass
'''


def gapfill_analysis(model_libsbml: Model, gapfill_params: dict[str: str], filename: str) -> Union[pd.DataFrame, tuple]:  # (Genbank) GFF file
    """Main function to infer gaps in a model by comparing the locus tags of the GeneProducts 
        to KEGG/BioCyc/both
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            
            gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
            
            filename (str): Path to output file for gapfill analysis result
            
        Returns:
            Case 'KEGG' - A data frame containing the columns 'bigg_id' 'locus_tag' 'EC' 'KEGG' 'name' 'GPR'
            
            Case 'BioCyc' - Five data frames:
            
                (1): Gap fill statistics with the columns 
                    'Missing entity' 'Total' 'Have BiGG ID' 'Can be added' 'Notes'
                    
                (2): Genes with the columns 
                    'locus_tag' 'protein_id' 'model_id' 'name'
                
                (3): Metabolites with the columns 
                    'bigg_id' 'name' 'BioCyc' 'compartment' 'Chemical Formula' 'InChI-Key' 'ChEBI' 'charge'
                    
                (4): Metabolites without BiGG ID with the columns 
                    'BioCyc' 'Chemical Formula' 'InChI-Key' 'ChEBI' ('charge')
                    
                (5): Reactions with the columns 
                    'bigg_id' 'name' 'BioCyc' 'locus_tag' 'Reactants' 'Products' 'EC' 'Fluxes' 'Spontaneous?' 
                    'bigg_reaction'
                    
            Case 'KEGG+BioCyc': Same output as cases 'KEGG' & 'BioCyc' 
                -> Data frame reactions contains additionally column 'KEGG'
    """
    colorama_init(autoreset=True)
    db_to_compare = gapfill_params['db_to_compare']
    result = None
    
    if db_to_compare not in ['KEGG', 'BioCyc', 'KEGG+BioCyc']:  # 'GFF', 
        print(f'{Fore.RED}To use the module gapfill the parameter of db_to_compare has to be set to one of the following' 
              + ' options:\n- \'KEGG\'\n- \'BioCyc\'\n- \'KEGG+BioCyc\'\nAdditionally, the required parameters' 
              + ' for each option need to be specified.\n- \'biggreactions\' and \'gapfill\' are required for all options.'
              + '\n- \'organismid\' is required only for the options \'KEGG\' and \'KEGG+BioCyc\'.\n- \'biocyc_tables\'' 
              + ' is only required for the options \'BioCyc\' and \'KEGG+BioCyc\'.')  # \n- \'GFF\'
        return
   
    if db_to_compare == 'KEGG':
        if gapfill_params['organismid']:
            missing_kegg = rga_kegg.kegg_gene_comp(model_libsbml, 
                                                   gapfill_params['organismid'], 
                                                   gapfill_params['gff_file']
                                                   )
            result = missing_kegg
        else:
            print(f'{Fore.RED}To use the KEGG comparison the specification of the organismid is obligatory.\n' +
                  'If there is no organismid available for your organism in KEGG, use one of the options \'BioCyc\' or \'GFF\'.')
        
    elif db_to_compare == 'BioCyc':
        missing_biocyc = rga_biocyc.biocyc_gene_comp(model_libsbml, 
                                                     gapfill_params['biocyc_files']
                                                     )
        result = missing_biocyc
        
        ''' Implement here call of function that can be used with lab strain/organism which is in no database contained
        elif db_to_compare == 'GFF':
            gff_genes = gff_gene_comp(model_libsbml, 
                                    gapfill_params['gff_file']
                                    )
            result = gff_genes
        '''
        
    elif db_to_compare == 'KEGG+BioCyc':
        missing_kegg_reacs = rga_kegg.kegg_gene_comp(model_libsbml, 
                                               gapfill_params['organismid'], 
                                               gapfill_params['gff_file']
                                               )
        missing_kegg_reacs.drop(['name', 'locus_tag', 'EC'], axis=1, inplace=True)
        missing_biocyc = rga_biocyc.biocyc_gene_comp(model_libsbml, 
                                                     gapfill_params['biocyc_files']
                                                     )
        stats, missing_biocyc_genes, missing_biocyc_metabs, missing_metabs_wo_BiGG_df, missing_biocyc_reacs = missing_biocyc
        missing_combined_reacs = missing_biocyc_reacs.merge(missing_kegg_reacs[['bigg_id', 'KEGG']], how='left', on='bigg_id')
        result = (stats, missing_biocyc_genes, missing_biocyc_metabs, missing_metabs_wo_BiGG_df, missing_combined_reacs, missing_kegg_reacs)
        
        
    if type(result) == tuple:
        with pd.ExcelWriter(filename) as writer:
            result[0].to_excel(writer, sheet_name='gap fill statistics', index=False)
            result[1].to_excel(writer, sheet_name='genes', index=False)
            result[2].to_excel(writer, sheet_name='metabolites', index=False)
            result[3].to_excel(writer, sheet_name='metabolites without BiGG IDs', index=False)
            result[4].to_excel(writer, sheet_name='reactions', index=False)
            if len(result) == 6:
                result[5].to_excel(writer, sheet_name='KEGG reactions', index=False)
    else:
        with pd.ExcelWriter(filename) as writer:
            result.to_excel(writer, sheet_name='KEGG reactions', index=False)
        
    return result
    
    
def gapfill_model(model_libsbml: Model, gapfill_analysis_result: Union[str, tuple]):
    """Main function to fill gaps in a model from a table
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            
            gapfill_analysis_result (str | tuple): Path to Excel file from gapfill_analysis | Tuple of pandas data frames obtained from gapfill_analysis
            
        Return:
            Gap filled model
    """
    model = model_libsbml
    missing_genes_df, missing_metabs_df, missing_reacs_df = None, None, None
    if type(gapfill_analysis_result) == tuple:  # Tuple of pandas dataframes from gapfill_analysis
        missing_genes_df = gapfill_analysis_result[1]
        missing_metabs_df = gapfill_analysis_result[2]
        missing_reacs_df = gapfill_analysis_result[4]
    else:  # Excel file from user-input
        with pd.ExcelFile(gapfill_analysis_result) as reader:
            gp_analysis_res = pd.read_excel(reader, sheet_name=['genes', 'metabolites', 'reactions'])
        missing_genes_df = gp_analysis_res.get('genes')
        missing_metabs_df = gp_analysis_res.get('metabolites')
        missing_reacs_df = gp_analysis_res.get('reactions')
        
    missing_genes_df = missing_genes_df.replace(np.nan, None)
    missing_metabs_df = missing_metabs_df.replace(np.nan, None)
    missing_reacs_df = missing_reacs_df.replace(np.nan, None)
    
    # (1) Add all missing genes needed for the missing reactions
    for _, row in missing_genes_df.iterrows():
        gp, model = create_gp(model_libsbml, row['model_id'], row['name'], row['locus_tag'], row['protein_id'])
    
    # (2) Add all missing metabolites needed for the missing reactions
    for _, row in missing_metabs_df.iterrows():
        sp, model = create_species(model_libsbml, row['bigg_id'], row['name'], row['compartment'], row['charge'], row['Chemical Formula'])
        if 'BioCyc' in missing_metabs_df.columns:
            biocyc_row = ast.literal_eval(str(row['BioCyc']))
            if biocyc_row:
                for biocyc_id in biocyc_row:
                    add_cv_term_metabolites(biocyc_id, 'BioCyc', sp)
                    add_cv_term_metabolites(biocyc_id, 'METACYC', sp)
                
        if 'InChI-Key' in missing_metabs_df.columns:
            if row['InChI-Key']:
                add_cv_term_metabolites(str(row['InChI-Key']), 'InChI-Key', sp)
                
        if 'ChEBI' in missing_metabs_df.columns:
            if row['ChEBI']:
                add_cv_term_metabolites(str(row['ChEBI']), 'ChEBI', sp)
    
    # (3) Add all missing reactions
    for _, row in missing_reacs_df.iterrows():
        reaction_dict = ast.literal_eval(str(row['bigg_reaction']))
        reactants = reaction_dict.get('reactants')
        products = reaction_dict.get('products')
        genes = ast.literal_eval(str(row['gene_product'])) if row['Spontaneous?'] != 'T' else 'G_spontaneous'
        compartment = row['compartment']
        compartment = compartment if compartment != 'exchange' else None
        reac, model = create_reaction(
            model=model_libsbml, reaction_id=row['bigg_id'], name=row['name'], reactants=reactants, 
            products=products, fluxes=ast.literal_eval(str(row['fluxes'])), compartment=compartment, genes=genes
            )
        if 'bigg_aliases' in missing_reacs_df.columns:
            bigg_aliases_row = ast.literal_eval(str(row['bigg_aliases']))
            if bigg_aliases_row:
                for bigg_id in bigg_aliases_row:
                    if bigg_id != reac.getId():
                        add_cv_term_reactions(bigg_id, 'BIGG', reac)
        
        if 'KEGG' in missing_reacs_df.columns:
            kegg_row = ast.literal_eval(str(row['KEGG']).replace('nan', 'None'))
            if kegg_row:
                for kegg_id in kegg_row:
                    add_cv_term_reactions(kegg_id, 'KEGG', reac)
                
        if 'BioCyc' in missing_reacs_df.columns:
            biocyc_row = ast.literal_eval(str(row['BioCyc']))
            if biocyc_row:
                for biocyc_id in biocyc_row:
                    add_cv_term_reactions(biocyc_id, 'BioCyc', reac)
                    add_cv_term_reactions(biocyc_id, 'METACYC', reac)
                    
        if 'EC' in missing_reacs_df.columns:
            ec_row = ast.literal_eval(str(row['EC']))
            if ec_row:
                for ec_num in ec_row:
                    add_cv_term_reactions(ec_num, 'EC', reac)
                         
    return model


def gapfill(
    model_libsbml: Model, gapfill_params: dict[str: str], filename: str
    ) -> Union[tuple[pd.DataFrame, Model], tuple[tuple, Model]]:
    """Main function to fill gaps in a model by comparing the locus tags of the GeneProducts to 
        KEGG/BioCyc/(Genbank) GFF file
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            
            gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
            
            filename (str): Path to output file for gapfill analysis result
            
            gapfill_model_out (str): Path where gapfilled model should be written to
            
        Return:
            Result from gapfill_analysis() and gap filled model
    """
    gapfill_analysis_result = gapfill_analysis(model_libsbml, gapfill_params, filename)
    model = gapfill_model(model_libsbml, gapfill_analysis_result)
    
    return gapfill_analysis_result, model
