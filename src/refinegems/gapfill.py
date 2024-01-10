#!/usr/bin/env python
"""The `gapfill` module can be used either with KEGG were you only need the KEGG organism ID or with BioCyc or with both (Options: 'KEGG', 'BioCyc', 'KEGG+BioCyc').
    For how to obtain the BioCyc tables look into the documentation under 'Filling gaps with refineGEMs' > 'Automated gap filling'.
    
    Run times: 
    
        * 'KEGG': ~ 2h
        * 'BioCyc': ~ 45mins - 1h
        * 'KEGG+BioCyc': ~ 3 - 4h  
"""
import ast
import math
from libsbml import Model as libModel
import refinegems.analysis_kegg as rga_kegg
import refinegems.analysis_biocyc as rga_biocyc
from refinegems.curate import update_annotations_from_others
from refinegems.cvterms import add_cv_term_metabolites, add_cv_term_reactions 
from refinegems.entities import create_gp, create_species, create_reaction
import pandas as pd
import numpy as np
from typing import Union
from colorama import init as colorama_init
from colorama import Fore

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel"

'''Skeleton for functions that could be used for a lab strain/organism which is in no database contained
def get_genes_from_gff():
    pass


def get_related_metabs_reactions_blast():
    pass


def gff_gene_comp():
    pass
'''


def gap_analysis(model_libsbml: libModel, gapfill_params: dict[str: str], filename: str) -> Union[pd.DataFrame, tuple]:  # (Genbank) GFF file
    """| Main function to infer gaps in a model by comparing the locus tags of the GeneProducts 
       | to KEGG/BioCyc/both

    Args:
        - model_libsbml (libModel): Model loaded with libSBML
        - gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings 
        - filename (str): Path to output file for gapfill analysis result
        
    Returns:
        - Case 'KEGG'
            pd.DataFrame: Table containing the columns 'bigg_id' 'locus_tag' 'EC' 'KEGG' 'name' 'GPR'
        - Case 'BioCyc'
            tuple: Five tables (1) - (4)
                (1) pd.DataFrame: Gap fill statistics with the columns 
                                    'Missing entity' 'Total' 'Have BiGG ID' 'Can be added' 'Notes'
                (2) pd.DataFrame: Genes with the columns 
                                    'locus_tag' 'protein_id' 'model_id' 'name'
                (3) pd.DataFrame: Metabolites with the columns 
                                    'bigg_id' 'name' 'BioCyc' 'compartment' 'Chemical Formula' 'InChI-Key' 'ChEBI' 'charge'  
                (4) pd.DataFrame: Reactions with the columns 
                                    'bigg_id' 'name' 'BioCyc' 'locus_tag' 'Reactants' 'Products' 'EC' 'Fluxes' 'Spontaneous?' 
                                    'bigg_reaction'
                
        - Case 'KEGG+BioCyc': 
            tuple: Five tables (1)-(4) from output of 'BioCyc' & (5) from output of 'KEGG'
                    -> Table reactions contains additionally column 'KEGG'
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
                  'If there is no organismid available for your organism in KEGG but an entry for your organism exists in BioCyc, use the option \'BioCyc\'.\n' +
                  'If no entry for your organism exists in KEGG and/or BioCyc, the gap analysis cannot be done.')
            # use one of the options \'BioCyc\' or \'GFF\'
        
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
        stats, missing_biocyc_genes, missing_biocyc_metabs, missing_biocyc_reacs = missing_biocyc 
        missing_combined_reacs = missing_biocyc_reacs.merge(missing_kegg_reacs[['bigg_id', 'KEGG']], how='left', on='bigg_id')
        result = (stats, missing_biocyc_genes, missing_biocyc_metabs, missing_combined_reacs, missing_kegg_reacs)
        
        
    if type(result) == tuple:
        with pd.ExcelWriter(f'{filename}.xlsx') as writer:
            result[0].to_excel(writer, sheet_name='gap fill statistics', index=False)
            result[1].to_excel(writer, sheet_name='genes', index=False)
            result[2].to_excel(writer, sheet_name='metabolites', index=False)
            result[3].to_excel(writer, sheet_name='reactions', index=False)
            if len(result) == 5:
                result[4].to_excel(writer, sheet_name='KEGG reactions', index=False)
    else:
        with pd.ExcelWriter(f'{filename}.xlsx') as writer:
            result.to_excel(writer, sheet_name='KEGG reactions', index=False)
        
    return result
    
    
def gapfill_model(model_libsbml: libModel, gap_analysis_result: Union[str, tuple]) -> libModel:
    """Main function to fill gaps in a model from a table

    Args:
        - model_libsbml (libModel): Model loaded with libSBML
        - gap_analysis_result (str|tuple): Path to Excel file from gap_analysis|Tuple of pd.DataFrames obtained from gap_analysis
        
    Returns:
        libModel: Gap filled model
    """
    model = model_libsbml
    missing_genes_df, missing_metabs_df, missing_reacs_df = None, None, None
    if type(gap_analysis_result) == tuple:  # Tuple of pandas dataframes from gap_analysis
        missing_genes_df = gap_analysis_result[1]
        missing_metabs_df = gap_analysis_result[2]
        missing_reacs_df = gap_analysis_result[3]
    else:  # Excel file from user-input
        with pd.ExcelFile(gap_analysis_result) as reader:
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
            biocyc_row = ast.literal_eval(str(row['BioCyc']).replace('nan', 'None'))
            if biocyc_row:
                for biocyc_id in biocyc_row:
                    if biocyc_id:
                        add_cv_term_metabolites(biocyc_id, 'BioCyc', sp)
                        add_cv_term_metabolites(biocyc_id, 'METACYC', sp)
                
        if 'InChI-Key' in missing_metabs_df.columns:
            inchi_key = str(row['InChI-Key'])
            inchi_key = inchi_key.removeprefix('InChIKey=') if 'InChIKey=' in inchi_key else inchi_key
            inchi_key = inchi_key if inchi_key != 'nan' else None  
            if inchi_key:
                add_cv_term_metabolites(inchi_key, 'InChI-Key', sp)

        if 'ChEBI' in missing_metabs_df.columns:
            chebi_value = row.get('ChEBI')
            if chebi_value is not None:
                chebi_id = str(int(chebi_value))
                add_cv_term_metabolites(chebi_id, 'ChEBI', sp)

    model = update_annotations_from_others(model)

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
                    if kegg_id:
                        add_cv_term_reactions(kegg_id, 'KEGG', reac)
                
        if 'BioCyc' in missing_reacs_df.columns:
            biocyc_row = ast.literal_eval(str(row['BioCyc']))
            if biocyc_row:
                for biocyc_id in biocyc_row:
                    add_cv_term_reactions(biocyc_id, 'BioCyc', reac)
                    add_cv_term_reactions(biocyc_id, 'METACYC', reac)
                    
        if 'EC' in missing_reacs_df.columns:
            ec_row = ast.literal_eval(str(row['EC']).replace('nan', 'None'))
            if ec_row:
                for ec_num in ec_row:
                    if ec_num:
                        add_cv_term_reactions(ec_num, 'EC', reac)
                         
    return model


def gapfill(
    model_libsbml: libModel, gapfill_params: dict[str: str], filename: str
    ) -> Union[tuple[pd.DataFrame, libModel], tuple[tuple, libModel]]:
    """| Main function to fill gaps in a model by comparing the locus tags of the GeneProducts to 
       | KEGG/BioCyc/(Genbank) GFF file

    Args:
        - model_libsbml (libModel): Model loaded with libSBML
        - gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
        - filename (str): Path to output file for gapfill analysis result
        - gapfill_model_out (str): Path where gapfilled model should be written to
        
    Returns:
        tuple: ``gap_analysis()`` table(s) (1) & libSBML model (2)
            (1) pd.DataFrame|tuple(pd.DataFrame): Result from function ``gap_analysis()``
            (2) libModel: Gap filled model
    """
    gap_analysis_result = gap_analysis(model_libsbml, gapfill_params, filename)
    model = gapfill_model(model_libsbml, gap_analysis_result)
    
    return gap_analysis_result, model
