#!/usr/bin/env python
from libsbml import *
import refinegems.analysis_db as rga
import refinegems.analysis_kegg as rga_kegg
import refinegems.analysis_biocyc as rga_biocyc
import refinegems.entities as rge
from colorama import init as colorama_init
from colorama import Fore

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


def get_genes_from_gff():
    pass


def get_related_metabs_reactions_blast():
    pass


def gff_gene_comp():
    pass


def gapfill_analysis(model_libsbml: Model, gapfill_params: dict[str: str]):
    """ Main function to gapfill a model with comparison to KEGG/BioCyc/(Genbank) GFF file
   
        Args:
            model_libsbml (Model): model loaded with libSBML
            gapfill_params (dict): Dictionary obtained from YAML file containing the parameter mappings
    """
    colorama_init(autoreset=True)
    db_to_compare = gapfill_params['db_to_compare']
    
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
                                                   gapfill_params['bigg_dbs'][0], 
                                                   gapfill_params['gff_file']
                                                   )
            return missing_kegg
        else:
            print(f'{Fore.RED}To use the KEGG comparison the specification of the organismid is obligatory.\n' +
                  'If there is no organismid available for your organism in KEGG, use one of the options \'BioCyc\' or \'GFF\'.')
        
    elif db_to_compare == 'BioCyc':
        missing_biocyc = rga_biocyc.biocyc_gene_comp(model_libsbml, 
                                                     gapfill_params['biocyc_files'],
                                                     gapfill_params['bigg_dbs']
                                                     )
        return missing_biocyc
        
    elif db_to_compare == 'GFF':
        gff_genes = gff_gene_comp(model_libsbml, 
                                  gapfill_params['bigg_dbs'], 
                                  gapfill_params['gff_file']
                                  )
        return gff_genes
        
    elif db_to_compare == 'KEGG+BioCyc':
        missing_kegg_reacs = rga_kegg.kegg_gene_comp(model_libsbml, 
                                               gapfill_params['organismid'], 
                                               gapfill_params['bigg_dbs'][0], 
                                               gapfill_params['gff_file']
                                               )
        missing_kegg_reacs.drop(['name', 'locus_tag', 'EC'], axis=1, inplace=True)
        missing_biocyc = rga_biocyc.biocyc_gene_comp(model_libsbml, 
                                                     gapfill_params['biocyc_files'],
                                                     gapfill_params['bigg_dbs']
                                                     )
        stats, missing_biocyc_genes, missing_biocyc_metabs, missing_metabs_wo_BiGG_df, missing_biocyc_reacs = missing_biocyc
        missing_combined_reacs = missing_biocyc_reacs.merge(missing_kegg_reacs, how='left', on='bigg_id')
        return (stats, missing_biocyc_genes, missing_biocyc_metabs, missing_metabs_wo_BiGG_df, missing_combined_reacs, missing_kegg_reacs)
