#!/usr/bin/env python
import re
import requests
import sqlite3
import sqlalchemy
import pandas as pd
pd.options.mode.chained_assignment = None # suppresses the pandas SettingWithCopyWarning; comment out before developing!!
import numpy as np
from refinegems.io import load_a_table_from_database
from refinegems.databases import PATH_TO_DB
from typing import Literal
from tqdm import tqdm
from ratelimit import limits, sleep_and_retry
from multiprocessing import Pool

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel and Tobias Fehrenbach"


ALL_BIGG_COMPARTMENTS_ONE_LETTER = ('c', 'e', 'p', 'm', 'x', 'r', 'v', 'n', 'g', 'u', 'l', 'h', 'f', 's', 'i', 'w', 'y')
ALL_BIGG_COMPARTMENTS_TWO_LETTER = ('im', 'cx', 'um', 'cm', 'mm')
BIGG_REACTIONS_URL = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
BIGG_METABOLITES_URL = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/'

COMPARTMENTS = ('c', 'e', 'p')


def get_search_regex(other_db: Literal['KEGG', 'BioCyc', 'SEED'], metabolites: bool) -> str:
    """Retrieves the search regex for BioCyc/KEGG/SEED to be used in the BiGG mapping

    Args:
        - other_db (Literal): Specifies if the search regex should be for BioCyc/KEGG/SEED
        - metabolites (bool): Is required if one wants to search for KEGG/SEED Compound IDs in the bigg_models_metabolites.txt
            
    Returns:
        str: Search regex
    """
    if other_db == 'BioCyc':
        return 'BioCyc: http://identifiers.org/biocyc/META:(.*?);'
    elif other_db == 'KEGG' or other_db == 'SEED':
        if metabolites:
            return f'{other_db} Compound: http://identifiers.org/{other_db.lower()}.compound/(.*?);'
        else:
            return f'{other_db} Reaction: http://identifiers.org/{other_db.lower()}.reaction/(.*?);'
        
        
def compare_ids(id1: str, id2: str) -> bool:
    """Compares two strings/IDs & Returns True if one string matches most of the other

    Args:
        - id1 (str): ID 1
        - id2 (str): ID 2

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
 

@sleep_and_retry
@limits(calls=10, period=1)
def get_reaction_compartment(bigg_id: str) -> str:
    """Retrieves the compatment(s) a reaction is hapening in from the BiGG reaction identifier
        via the metabolites

    Args:
        - bigg_id (str): BiGG reaction identifier

    Returns:
        (1) str: 
            - Compartment of the provided reaction if reaction in single compartment
            - 'exchange' if reaction in multiple compartments
        (2) np.nan: 'NaN' if one of the found compartments is not in COMPARTMENTS
    """
    
    metabs_from_reac = requests.get(BIGG_REACTIONS_URL + bigg_id, allow_redirects=False).json()['metabolites']
            
    comps = [comp_dict.get('compartment_bigg_id') for comp_dict in metabs_from_reac]  # Get all compartments for reaction
    contained_in_compartments = [(comp in COMPARTMENTS) for comp in comps]  # Get True for correct compartment        
    if not all(contained_in_compartments):  # At least one compartment not correct
        return np.nan
    else:  # All compartments correct
        if len(set(comps)) == 1:  # Set of found compartments of reaction = 1: Reaction happens in one compartment
            return comps[0]
        else:  # Not so important but do not remove reaction as reaction in correct compartments
            return 'exchange'  # Probably exchange reaction


def keep_only_reactions_in_certain_compartments(complete_df: pd.DataFrame) -> pd.DataFrame:
    """Extracts all possible BiGG ID variations from database for a BiGG reaction ID, gets the metabolite compartments
        & returns table containing only reactions which happen in one of the provided compartments
        
    Args:
        - complete_df (pd.DataFrame): Table containing at least the columns 'bigg_id' & 'KEGG'/'BioCyc'
        
    Returns:
        pd.DataFrame: Table containing reactions & their compartments
    """
    tqdm.pandas()
    db = 'KEGG' if 'KEGG' in complete_df.columns else 'BioCyc'
    complete_df = complete_df[['bigg_id', db]]  # Remove all unnecessary columns
    
    # (1) Find all occurrencs of a BiGG reaction ID in bigg_reactions table in database
    def get_all_similar_bigg_ids(bigg_id_in: str) -> list[str]:
        
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
    # get_react_compartment moved to outer scope due to multiprocessing Pool (see line 95)
    def multi_get_reaction_compartment(complete_df: pd.DataFrame) -> list:
        """Takes a dataframe and runs get_reaction_compartment() in multiple
            processes on the 'bigg_id' column

        Args:
            - complete_df (pd.DataFrame): Table containing at least the columns 'bigg_id' & 'KEGG'/'BioCyc'

        Returns:
            list: List of compartments
        """
        with Pool() as pool:
            results = []
            for out in tqdm(pool.imap(get_reaction_compartment, complete_df.loc[:, "bigg_id"], chunksize=20), total=len(complete_df)):
                results.append(out)

        return results

    # Connect to database & get similar IDs (1)
    print('Getting all similar IDs...')
    con = sqlite3.connect(PATH_TO_DB)  # Open connection to database
    complete_df.loc[:, 'bigg_id_list'] = complete_df.loc[:, 'bigg_id'].progress_map(get_all_similar_bigg_ids)
    # complete_df.progress_apply(get_all_similar_bigg_ids, axis=1)
    con.close()  # Close connection to database

    # Adjust table to contain one BiGG ID per row from bigg_id_list (1)
    complete_df.loc[:, 'id_group'] = complete_df['bigg_id'].ne(complete_df['bigg_id'].shift()).cumsum()  # Group similar IDs
    complete_df.drop(labels='bigg_id', axis=1, inplace=True)  # Drop 'bigg_id' as no longer required
    complete_df = complete_df.explode('bigg_id_list', ignore_index=True)  # Expand 'bigg_id_list' column
    complete_df.rename(columns={'bigg_id_list': 'bigg_id'}, inplace=True)  # Rename 'bigg_id_list' to 'bigg_id'

    # (2) Get all compartments for each reaction from BiGG database API
    print(f'Getting all IDs with correct compartment {COMPARTMENTS}...')
    results = multi_get_reaction_compartment(complete_df)
    complete_df["compartment"] = results

    # complete_df.progress_apply(get_reaction_compartment, axis=1)  # (2)

    # (3) Remove reactions with compartment = NaN
    complete_df.dropna(subset=['compartment'], inplace=True)

    return complete_df


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_bigg2other_db(other_db: Literal['KEGG', 'BioCyc', 'SEED'], metabolites: bool=False) -> pd.DataFrame:
    """Uses list of BiGG reactions/metabolites to get a mapping from BiGG to KEGG/BioCyc Id

    Args:
        - other_db (Literal): Set to 'KEGG'/'BioCyc'/'SEED' to map KEGG/BioCyc/SEED IDs to BiGG IDs
        - metabolites (bool): Set to True to map other_db IDs to BiGG IDs for metabolites

    Returns:
        pd.DataFrame: Table containing BiGG Ids with corresponding KEGG/BioCyc/SEED Ids
    """
    
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
    
    bigg_db_df[other_db] = bigg_db_df.apply(
        lambda row: find_other_db(row['database_links']), axis=1)
    bigg_db_df = bigg_db_df.explode(other_db, ignore_index=True)
    
    if not metabolites:
        bigg_db_df = keep_only_reactions_in_certain_compartments(bigg_db_df)
        
    bigg_df = bigg_db_df[['bigg_id', other_db]] if metabolites else bigg_db_df[['bigg_id', other_db, 'compartment', 'id_group']]

    return bigg_df
 
 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_bigg_model(complete_df: pd.DataFrame, model_entities: pd.DataFrame, metabolites: bool=False) -> pd.DataFrame:
    """Compares missing entities obtained through genes extracted via KEGG/BioCyc to entities in the model
        Needed to back check previous comparisons.

    Args:
        - complete_df (pd.DataFrame): Table that contains KEGG/BioCyc Id, BiGG Id & more
        - model_entities (pd.DataFrame): BiGG Ids of entities in the model 
        - metabolites (bool): True if names of metabolites should be added, otherwise false

    Returns:
        pd.DataFrame: Table containing entities present in KEGG/BioCyc but not in the model
    """
    db = 'KEGG' if 'KEGG' in complete_df.columns else 'BioCyc'  # Find out which database was used
    
    # Get only IDs that are not in model
    mapp = complete_df.set_index('bigg_id')
    entities = model_entities.set_index('bigg_id')
    entities_missing_in_model = mapp[~mapp.index.isin(
        entities.index)].reset_index()
    
    db_ids = entities_missing_in_model.groupby('bigg_id')[db].agg(set)  # Get a set of all BioCyc/KEGG IDs belonging to one BiGG ID
    
    # Add set of BioCyc/KEGG IDs belonging to one BiGG ID to the dataframe
    entities_missing_in_model.set_index('bigg_id', inplace=True)
    entities_missing_in_model.loc[:, db] = db_ids
    entities_missing_in_model.reset_index(inplace=True)
    
    if 'id_group' in entities_missing_in_model.columns:  # Remove reaction ID duplicates but keep all related BiGG & BioCyc/KEGG IDs in a list
        aliases = entities_missing_in_model.groupby(['compartment', 'id_group'])['bigg_id'].agg(set)  # Get a set of the 'duplicated' BiGG reaction IDs -> aliases
        entities_missing_in_model.drop_duplicates(['compartment', 'id_group'], inplace=True, ignore_index=True)  # Drop duplicates where compartments & id_group same
        
        # Add set of BiGG ID aliases to the dataframe
        entities_missing_in_model.set_index(['compartment', 'id_group'], inplace=True)
        entities_missing_in_model.loc[:, 'bigg_aliases'] = aliases
        entities_missing_in_model.reset_index(inplace=True)
        
        entities_missing_in_model.drop(labels='id_group', axis=1, inplace=True)  # id_group is not longer necessary
        
    entities_missing_in_model.drop_duplicates(subset='bigg_id', inplace=True, ignore_index=True)  # Remove BiGG ID duplicates
    
    # Add name column to dataframe
    def get_name_from_bigg(bigg_id: str):
        bigg_db = 'bigg_metabolites' if metabolites else 'bigg_reactions'
        query = f"SELECT name FROM {bigg_db} WHERE bigg_id=\'{bigg_id}\'"
        name_from_bigg = con.execute(query).fetchone()[0]
        return name_from_bigg
    
    con = sqlite3.connect(PATH_TO_DB)  # Open connection to database
    entities_missing_in_model['name'] = entities_missing_in_model['bigg_id'].map(get_name_from_bigg)
    con.close()
    
    # Add compartment ID to all BiGG metabolites
    if metabolites:
        def get_compartment_from_id(bigg_id: str):
            compartment = bigg_id[-1]
            return compartment if compartment in COMPARTMENTS else np.nan  # To filter the incorrect compartments out
        
        entities_missing_in_model['compartment'] = entities_missing_in_model.apply(
            lambda row: get_compartment_from_id(row['bigg_id']), axis=1)
        entities_missing_in_model.dropna(subset=['compartment'], inplace=True)  # Drop all BiGG metabolite IDs which have no valid compartment
    
    return entities_missing_in_model


def add_stoichiometric_values_to_reacs(missing_reacs: pd.DataFrame) -> pd.DataFrame:
    """Adds for each reaction a dictionary containing the reactants & products as dictionaries with the BiGG Metabolite 
        ID as key & the respective absolute stoichiometric value as value
        
    Args:
        - missing_reacs (pd.DataFrame): Table containing missing reactions (Only requires a column containing BiGG IDs)
            
    Returns:
        pd.DataFrame: Table where for each BiGG reaction ID a dictionary containing reactants & products exists 
    """
    
    def get_reactants_and_products_dicts(reaction_id: str) -> list[dict]:
        reactants = {}
        products = {}
        
        metabs_from_reac = requests.get(BIGG_REACTIONS_URL + reaction_id).json()['metabolites']

        for compound_dict in metabs_from_reac:
            complete_bigg_id = None
            if compound_dict.get('compartment_bigg_id'):
                complete_bigg_id = f"{compound_dict.get('bigg_id')}_{compound_dict.get('compartment_bigg_id')}"
            else:
                complete_bigg_id = compound_dict.get('bigg_id')
            if compound_dict.get('stoichiometry') < 0:
                reactants[complete_bigg_id] = abs(compound_dict.get('stoichiometry'))
            elif compound_dict.get('stoichiometry') > 0:
                products[complete_bigg_id] = abs(compound_dict.get('stoichiometry'))
                
        return str({'reactants': reactants, 'products': products})
                
    missing_reacs['bigg_reaction']= missing_reacs.apply(
        lambda row: get_reactants_and_products_dicts(str(row['bigg_id'])), axis=1)  #, missing_reacs['bigg_products'], result_type='expand'
      
    return missing_reacs
 