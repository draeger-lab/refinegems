#!/usr/bin/env python
import re
import requests
import pandas as pd
from libsbml import *
from bioservices.kegg import KEGG
from typing import Literal

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


def get_search_regex(other_db: Literal['KEGG', 'BioCyc'], metabolites: bool) -> str:
    """Retrieves the search regex for BioCyc/KEGG to be used in the BiGG mapping

        Args:
            other_db (Literal): specifies if the search regex should be for BioCyc or KEGG
            metabolites (bool): is required if one wants to search for KEGG Compound IDs in the bigg_models_metabolites.txt
            
        Returns:
            str: search regex
    """
    if other_db == 'BioCyc':
        return 'BioCyc: http://identifiers.org/biocyc/META:(.*?);'
    elif other_db == 'KEGG':
        if metabolites:
            return 'KEGG Compound: http://identifiers.org/kegg.compound/(.*?);'
        else:
            return 'KEGG Reaction: http://identifiers.org/kegg.reaction/(.*?);'
    
 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_bigg2other_db(bigg_db: str, other_db: Literal['KEGG', 'BioCyc'], metabolites: bool=False) -> pd.DataFrame:
    """Uses list of BiGG reactions/metabolites to get a mapping from BiGG to KEGG/BioCyc Id

    Args:
        bigg_db (Str): path to file containing BiGG database
                       (Either: bigg_models_metabolites.txt (metabolites=True) 
                        Or: bigg_models_reactions.txt (metabolites=False))
        other_db (Literal): Set to 'KEGG'/'BioCyc' to map KEGG/BioCyc IDs to BiGG IDs
        metabolites (bool): set to True to map other_db IDs to BiGG IDs for metabolites

    Returns:
        df: table containing BiGG Ids with corresponding KEGG/BioCyc Ids
    """
    db_search_regex = get_search_regex(other_db, metabolites)
    # make the download of biggreactions/biggmetabolites possible to maintain database
    bigg_db_df = pd.read_csv(
        bigg_db,
        sep='\t').drop(
        'model_list',
        axis=1).dropna()

    def find_other_db(database_links: str):
        m = re.findall(
            db_search_regex,
            database_links)
        if m:
            return m
        else:
            return None
    
    def get_compartment_from_id(bigg_id: str):
        compartment = bigg_id[-1]
        return compartment if compartment in ['c', 'e', 'p'] else 'c'

    bigg_db_df[other_db] = bigg_db_df.apply(
        lambda row: find_other_db(row['database_links']), axis=1)
    if metabolites:
        bigg_db_df['compartment'] = bigg_db_df.apply(
            lambda row: get_compartment_from_id(row['bigg_id']), axis=1)
    bigg_db_df = bigg_db_df.explode(other_db, ignore_index=True)

    return bigg_db_df[['bigg_id', 'name', other_db, 'compartment']] if metabolites else bigg_db_df[['bigg_id', 'name', other_db]]
 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_bigg_model(complete_df: pd.DataFrame, model_entities: pd.DataFrame):
    """Compares missing entities obtained through genes extracted via KEGG/BioCyc to entities in the model
        Needed to back check previous comparisons.

    Args:
        complete_df (df): pandas dataframe that contains BioCyc Id, BiGG Id & more
        model_entities (df): BiGG Ids of entities in the model 

    Returns:
        df: table containing entities present in KEGG/BioCyc but not in the model
    """
    db = 'KEGG' if 'KEGG' in complete_df.columns else 'BioCyc'
    
    mapp = complete_df.set_index('bigg_id')
    entities = model_entities.set_index('bigg_id')

    entities_missing_in_model = mapp[~mapp.index.isin(
        entities.index)].reset_index()

    ambig_db_id = complete_df.loc[complete_df.duplicated(
        subset=[db], keep='first')]
    if 'EC' in ambig_db_id.columns:
        ambig_db_id = ambig_db_id.set_index(db).drop(
            ['locus_tag', 'EC'], axis=1).sort_index()
    else:
        ambig_db_id.set_index(db).sort_index(inplace=True)

    ambig = ambig_db_id.set_index('bigg_id')
    miss = entities_missing_in_model.set_index('bigg_id')

    entities_missing_in_model_non_dup = miss[~miss.index.isin(
        ambig.index)].reset_index()

    return entities_missing_in_model_non_dup


def add_stoichiometric_values_to_reacs(missing_reacs: pd.DataFrame) -> pd.DataFrame:
    
    def get_reactants_and_products_dicts(reaction_id: str) -> list[dict]:
        reactants = {}
        products = {}
        
        reac_bigg_url = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
        metabs_from_reac = requests.get(reac_bigg_url + reaction_id).json()['metabolites']

        for compound_dict in metabs_from_reac:
            if compound_dict.get('stoichiometry') < 0:
                reactants[compound_dict.get('bigg_id')] = abs(compound_dict.get('stoichiometry'))
            elif compound_dict.get('stoichiometry') > 0:
                products[compound_dict.get('bigg_id')] = abs(compound_dict.get('stoichiometry'))
                
        return str({'reactants': reactants, 'products': products})
                
    missing_reacs['bigg_reaction']= missing_reacs.apply(
        lambda row: get_reactants_and_products_dicts(str(row['bigg_id'])), axis=1)  #, missing_reacs['bigg_products'], result_type='expand'
      
    return missing_reacs  
 