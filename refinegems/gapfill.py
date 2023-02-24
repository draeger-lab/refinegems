#!/usr/bin/env python
from libsbml import *
import refinegems.analysis_kegg as rga_kegg
import refinegems.analysis_biocyc as rga_biocyc
import pandas as pd
from typing import Union
from colorama import init as colorama_init
from colorama import Fore

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


def get_genes_from_gff():
    pass


def get_related_metabs_reactions_blast():
    pass


def gff_gene_comp():
    pass


def gapfill_analysis(model_libsbml: Model, gapfill_params: dict[str: str], filename: str) -> Union[pd.DataFrame, tuple]:
    """Main function to infer gaps in a model by comparing the locus tags of the GeneProducts 
        to KEGG/BioCyc/(Genbank) GFF file
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
            filename (str): Path to output file for BioCyc
            
        Returns:
            Case 'KEGG' - A data frame containing the columns 'bigg_id' 'locus_tag' 'EC' 'KEGG' 'name' 'GPR'
            Case 'BioCyc' - Five data frames:
                (1): Gap fill statistics with the columns 
                    'Missing entity' 'Total' 'Have BiGG ID' 'Can be added' 'Notes'
                (2): Genes with the columns 'locus_tag' 'protein_id' 'model_id' 'name'
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
    
    if db_to_compare not in ['KEGG', 'BioCyc', 'GFF', 'KEGG+BioCyc']:
        print(f'{Fore.RED}To use the module gapfill the parameter of db_to_compare has to be set to one of the following' 
              + ' options:\n- \'KEGG\'\n- \'BioCyc\'\n- \'GFF\'\n- \'KEGG+BioCyc\'\nAdditionally, the required parameters' 
              + ' for each option need to be specified.\n- \'biggreactions\' and \'gapfill\' are required for all options.'
              + '\n- \'organismid\' is required only for the options \'KEGG\' and \'KEGG+BioCyc\'.\n- \'biocyc_tables\'' 
              + ' is only required for the options \'BioCyc\' and \'KEGG+BioCyc\'.')
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
        
    elif db_to_compare == 'GFF':
        gff_genes = gff_gene_comp(model_libsbml, 
                                  gapfill_params['gff_file']
                                  )
        result = gff_genes
        
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
        missing_combined_reacs = missing_biocyc_reacs.merge(missing_kegg_reacs, how='left', on='bigg_id')
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
            gapfill_analysis_result (str | tuple): Path to Excel file from gapfill_analysis | Tuple of pandas data 
                frames obtained from gapfill_analysis
    """
    if type(gapfill_analysis_result) == tuple:  # Tuple of pandas dataframes from gapfill_analysis
        pass  # Handle file
    else:  # Excel file from user-input
        pass  # Handle file
    pass


def gapfill(model_libsbml: Model, gapfill_params: dict[str: str]) -> Union[pd.DataFrame, tuple]:
    """Main function to fill gaps in a model by comparing the locus tags of the GeneProducts to 
        KEGG/BioCyc/(Genbank) GFF file
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
            
        Return:
            Result from gapfill_analysis()
    """
    gapfill_analysis_result = gapfill_analysis(model_libsbml, gapfill_params)
    gapfill_model(model_libsbml, gapfill_analysis_result)
    
    return gapfill_analysis_result
