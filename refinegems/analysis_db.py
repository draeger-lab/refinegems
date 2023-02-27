#!/usr/bin/env python
import re
import requests
import sqlite3
import pandas as pd
import numpy as np
from refinegems.io import load_a_table_from_database
from refinegems.databases import PATH_TO_DB
from libsbml import *
from typing import Literal
from tqdm import tqdm
from ratelimit import limits, sleep_and_retry


__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


ALL_BIGG_COMPARTMENTS_ONE_LETTER = ('c', 'e', 'p', 'm', 'x', 'r', 'v', 'n', 'g', 'u', 'l', 'h', 'f', 's', 'i', 'w', 'y')
ALL_BIGG_COMPARTMENTS_TWO_LETTER = ('im', 'cx', 'um', 'cm', 'mm')


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
        
        
def compare_ids(id1: str, id2: str) -> bool:
    """Compares two strings/IDs & Returns True if one string matches most of the other

    Args:
        id1 (str): ID 1
        id2 (str): ID 2

    Returns:
        bool: Indicates if most of one string contained in the other
    """
    id1_split, id2_split, id1_single_comp, id2_single_comp, id1_comp, id2_comp = None, None, None, None, None, None
    
    if '_' in id1: id1_split = re.split('_([a-zA-Z]|[0-9])$', id1)[0]
    if '_' in id2: id2_split = re.split('_([a-zA-Z]|[0-9])$', id2)[0]
    if id1.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): id1_single_comp = id1[:-1]
    if id2.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): id2_single_comp = id2[:-1]
    if id1.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): id1_comp = id1[:-2]
    if id2.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): id2_comp = id2[:-2]
    
    similar_ids = False
    if id1 == id2: similar_ids = True  # Both IDs are same
    
    elif id1_split and id2_split and (id1_split == id2_split): similar_ids = True # Both IDs are same but from different compartments
    elif id2_split and (id1 == id2_split): similar_ids = True # - "" -
    elif id1_split and (id1_split == id2): similar_ids = True # - "" -
    
    elif id1_single_comp and id2_single_comp and (id1_single_comp == id2_single_comp): similar_ids = True
    elif id2_single_comp and (id1 == id2_single_comp): similar_ids = True
    elif id1_single_comp and (id1_single_comp == id2_single_comp): similar_ids = True 
    
    elif id1_comp and id2_comp and (id1_comp == id2_comp): similar_ids = True
    elif id2_comp and (id1 == id2_comp): similar_ids = True
    elif id1_comp and (id1_comp == id2): similar_ids = True
    
    elif id1_split and id2_single_comp and (id1_split == id2_single_comp): similar_ids = True
    elif id2_split and id1_single_comp and (id1_single_comp == id2_split): similar_ids = True
    
    elif id1_split and id2_comp and (id1_split == id2_comp): similar_ids = True
    elif id2_split and id1_comp and (id1_comp == id2_split): similar_ids = True
    
    elif id1_comp and id2_single_comp and (id1_comp == id2_single_comp): similar_ids = True
    elif id2_comp and id1_single_comp and (id1_single_comp == id2_comp): similar_ids = True
    
    else: similar_ids = False

    return similar_ids

def keep_only_reactions_in_certain_compartments(complete_df: pd.DataFrame, compartments: tuple[str]) -> pd.DataFrame:
    """Extracts all possible BiGG ID variations from database for a BiGG reaction ID, gets the metabolite compartments
        & returns a dataframe containing only reactions which happen in one of the provided compartments
        
        Args:
            complete_df (DataFrame): A pandas dataframe containing at least the columns 'bigg_id' & 'KEGG'/'BioCyc'
            compartments (tuple): A tuple of BiGG compartment identifiers for which reaction IDs should be kept
        
        Returns:
            df: A pandas dataframe containing reactions & their compartments
    """
    tqdm.pandas()
    db = 'KEGG' if 'KEGG' in complete_df.columns else 'BioCyc'
    complete_df = complete_df[['bigg_id', 'name', db]]  # Remove all unnecessary columns
    
    # (1) Find all occurrencs of a BiGG reaction ID in bigg_reactions table in database
    def get_all_similar_bigg_ids(row: pd.Series) -> list[str]:
        bigg_id_in = str(row['bigg_id'])
        
        if '_' in bigg_id_in: bigg_id = re.split('_([a-zA-Z]|[0-9])$', bigg_id_in)[0]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): bigg_id = bigg_id_in[:-1]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): bigg_id = bigg_id_in[:-2]
        else: bigg_id = bigg_id_in
        
        query = f"SELECT bigg_id, INSTR(bigg_id, '{bigg_id}') bi FROM bigg_reactions WHERE bi > 0"
        result = con.execute(query).fetchall()
        result = [result_tuple[0] for result_tuple in result] if result else [bigg_id_in]
        result = [res for res in result if compare_ids(bigg_id, res)]
        return result
    
    # (2) Use list of all BiGG IDs obtained from database table bigg_reactions to get 'metabolites'
    @sleep_and_retry
    @limits(calls=10, period=1)
    def get_reaction_compartment(row: pd.Series) -> str:
        bigg_id = str(row['bigg_id'])
        
        reac_bigg_url = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
        metabs_from_reac = requests.get(reac_bigg_url + bigg_id, allow_redirects=False).json()['metabolites']
                
        comps = [comp_dict.get('compartment_bigg_id') for comp_dict in metabs_from_reac]  # Get all compartments for reaction
        contained_in_compartments = [True for comp in comps if comp in compartments]  # Get True for correct compartment        
        if not any(contained_in_compartments):  # At least one compartment not correct
            return np.nan
        else:  # All compartments correct
            if len(set(comps)) == 1:  # Set of found compartments of reaction = 1: Reaction happens in one compartment
                return comps[0]
            else:  # Not so important but do not remove reaction as reaction in correct compartments
                return 'exchange'  # Probably exchange reaction
    
    # Connect to database & get similar IDs (1)
    print('Getting all similar IDs...')
    con = sqlite3.connect(PATH_TO_DB)  # Open connection to database
    complete_df.loc[:,'bigg_id_list'] = complete_df.progress_apply(get_all_similar_bigg_ids, axis=1)
    con.close()  # Close connection to database
    
    # Adjust table to contain one BiGG ID per row from bigg_id_list (1)
    complete_df.loc[:, 'id_group'] = complete_df['bigg_id'].ne(complete_df['bigg_id'].shift()).cumsum()  # Group similar IDs
    complete_df.drop(labels='bigg_id', axis=1, inplace=True)  # Drop 'bigg_id' as no longer required
    complete_df = complete_df.explode('bigg_id_list', ignore_index=True)  # Expand 'bigg_id_list' column
    complete_df.rename(columns={'bigg_id_list': 'bigg_id'}, inplace=True)  # Rename 'bigg_id_list' to 'bigg_id'
    
    # (2) Get all compartments for each reaction from BiGG database API
    print(f'Getting all IDs with correct compartment {compartments}...')
    complete_df.loc[:, 'compartment'] = complete_df.progress_apply(get_reaction_compartment, axis=1)  # (2)
    
    # (3) Remove reactions with compartment = NaN
    complete_df.dropna(subset=['compartment'], inplace=True)
        
    return complete_df

 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_bigg2other_db(other_db: Literal['KEGG', 'BioCyc'], metabolites: bool=False) -> pd.DataFrame:
    """Uses list of BiGG reactions/metabolites to get a mapping from BiGG to KEGG/BioCyc Id

    Args:
        other_db (Literal): Set to 'KEGG'/'BioCyc' to map KEGG/BioCyc IDs to BiGG IDs
        metabolites (bool): set to True to map other_db IDs to BiGG IDs for metabolites

    Returns:
        df: table containing BiGG Ids with corresponding KEGG/BioCyc Ids
    """
    compartments = ('c', 'e', 'p')
    
    # Get only rows with BioCyc/KEGG entries
    db_table_name = 'bigg_metabolites' if metabolites else 'bigg_reactions'
    reaction_or_compound = 'Compound' if metabolites else 'Reaction'
    other_db_query = other_db if other_db == 'BioCyc' else ' '.join([other_db, reaction_or_compound])
    bigg_db_query = f"SELECT *, INSTR(database_links, '{other_db_query}:') o_db FROM {db_table_name} WHERE o_db > 0"
    bigg_db_df = load_a_table_from_database(bigg_db_query)
    
    db_search_regex = get_search_regex(other_db, metabolites)
    
    def find_other_db(database_links: str):
        m = re.findall(
            db_search_regex,
            str(database_links))
        if m:
            return m
        else:
            return None
    
    def get_compartment_from_id(bigg_id: str):
        compartment = bigg_id[-1]
        return compartment if compartment in compartments else np.nan  # To filter the incorrect compartments out
    
    bigg_db_df[other_db] = bigg_db_df.apply(
        lambda row: find_other_db(row['database_links']), axis=1)
    bigg_db_df = bigg_db_df.explode(other_db, ignore_index=True)
    if metabolites:
        bigg_db_df['compartment'] = bigg_db_df.apply(
            lambda row: get_compartment_from_id(row['bigg_id']), axis=1)
        bigg_db_df.dropna(subset=['compartment'], inplace=True)  # Drop all BiGG metabolite IDs which have no valid compartment
    else:
        bigg_db_df = keep_only_reactions_in_certain_compartments(bigg_db_df, compartments)
        
    bigg_df = bigg_db_df[['bigg_id', 'name', other_db, 'compartment']] if metabolites else bigg_db_df[['bigg_id', 'name', other_db, 'compartment', 'id_group']]

    return bigg_df
 
 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_bigg_model(complete_df: pd.DataFrame, model_entities: pd.DataFrame) -> pd.DataFrame:
    """Compares missing entities obtained through genes extracted via KEGG/BioCyc to entities in the model
        Needed to back check previous comparisons.

    Args:
        complete_df (df): pandas dataframe that contains BioCyc Id, BiGG Id & more
        model_entities (df): BiGG Ids of entities in the model 

    Returns:
        df: table containing entities present in KEGG/BioCyc but not in the model
    """
    db = 'KEGG' if 'KEGG' in complete_df.columns else 'BioCyc'  # Find out which database was used
    
    # Get only IDs that are not in model
    mapp = complete_df.set_index('bigg_id')
    entities = model_entities.set_index('bigg_id')
    entities_missing_in_model = mapp[~mapp.index.isin(
        entities.index)].reset_index()
    
    if 'id_group' in entities_missing_in_model.columns:  # Remove reaction ID duplicates but keep all realted BiGG & BioCyc IDs in a list
        aliases = entities_missing_in_model.groupby(['compartment', 'id_group'])['bigg_id'].agg(tuple)  # Get a list of the 'duplicated' BiGG reaction IDs -> aliases
        db_ids = entities_missing_in_model.groupby(['compartment', 'id_group'])[db].agg(tuple)  # Get a list of all BioCyc/KEGG IDs belonging to one BiGG reaction ID
        entities_missing_in_model.drop_duplicates(['compartment', 'id_group'], inplace=True)  # Drop duplicates where compartments & id_group same
        # Add lists to the dataframe
        entities_missing_in_model.set_index(['compartment', 'id_group'], inplace=True)
        entities_missing_in_model.loc[:, 'bigg_aliases'] = aliases
        entities_missing_in_model.loc[:, db] = db_ids
        entities_missing_in_model.reset_index(inplace=True)
        entities_missing_in_model.drop(labels='id_group', axis=1, inplace=True)
        
    entities_missing_in_model.drop_duplicates(subset=['bigg_id', db], inplace=True, ignore_index=True)  # Remove ID duplicates
    return entities_missing_in_model


def add_stoichiometric_values_to_reacs(missing_reacs: pd.DataFrame) -> pd.DataFrame:
    """Adds for each reaction a dictionary containing the reactants & products as dictionaries with the BiGG Metabolite 
        ID as key & the respective absolute stoichiometric value as value
        
        Args:
            missing_reacs (df): A dataframe containing missing reactions (Only requires a column containing BiGG IDs)
            
        Returns:
            df: A table where for each BiGG reaction ID a dictionary containing reactants & products exists 
    """
    
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
 