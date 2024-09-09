#!/usr/bin/env python
"""Access information from different databases or compare a model or model entities
with them. This module provides variables and function for accessing databases 
for better model curation and annotation.

The following databases have functionalities implemented:

- BiGG
- ChEBI
- KEGG
- ModelSEED
- NCBI
- UniProt
"""

__author__ = """Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune, 
            Jan-Philipp Leusch and Tobias Fehrenbach"""

############################################################################
# requirements
############################################################################

import cobra
import io
import libchebipy
import math
import numpy as np
import pandas as pd
import re
import requests
import sqlite3
import xmltodict

from Bio import Entrez, SeqIO
from Bio.KEGG import REST, Gene, Enzyme
from bioservices.kegg import KEGG
from cobra import Model as cobraModel
from multiprocessing import Pool
from ratelimit import limits, sleep_and_retry
from typing import Literal
from tqdm import tqdm

tqdm.pandas()
pd.options.mode.chained_assignment = None # suppresses the pandas SettingWithCopyWarning; comment out before developing!!

from .connections import run_DIAMOND_blastp, filter_DIAMOND_blastp_results
from .databases import PATH_TO_DB
from .io import load_a_table_from_database, create_missing_genes_protein_fasta
from .util import VALID_COMPARTMENTS

############################################################################
# variables
############################################################################

# database urls
# -------------

BIGG_REACTIONS_URL = 'http://bigg.ucsd.edu/api/v2/universal/reactions/' #: :meta: 
BIGG_METABOLITES_URL = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/' #: :meta: 

# Compartments in BiGG namespace
# ------------------------------

ALL_BIGG_COMPARTMENTS_ONE_LETTER = ('c', 'e', 'p', 'm', 'x', 'r', 'v', 'n', 'g', 'u', 'l', 'h', 'f', 's', 'i', 'w', 'y') #: :meta: 
ALL_BIGG_COMPARTMENTS_TWO_LETTER = ('im', 'cx', 'um', 'cm', 'mm') #: :meta: 

############################################################################
# functions
############################################################################

# BiGG
# ----
"""Map an ID to information in BiGG or compare model entities to BiGG.
"""

# @NOTE: Add to reaction handling for `build_reac_bigg`?
def add_stoichiometric_values_to_reacs_from_bigg(missing_reacs: pd.DataFrame) -> pd.DataFrame:
    """Adds for each reaction a dictionary containing the reactants & products as dictionaries with the BiGG Metabolite 
    ID as key & the respective absolute stoichiometric value as value
        
    Args:
        - missing_reacs (pd.DataFrame): 
            Table containing missing reactions (Only requires a column containing BiGG IDs)
            
    Returns:
        pd.DataFrame: 
            Table where for each BiGG reaction ID a dictionary containing reactants & products exists 
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


# @NOTE: Add to reaction handling for `build_reac_bigg`?
@sleep_and_retry
@limits(calls=10, period=1)
def get_reaction_compartment_from_bigg(bigg_id: str) -> str:
    """Retrieves the compatment(s) a reaction is hapening in from the BiGG reaction identifier via the metabolites

    Args:
        - bigg_id (str): 
            BiGG reaction identifier

    Returns:

        (1) Case ``compartment in VALID_COMPARTMENTS.keys()``
                str: 
                    Either 
                    
                    - Compartment of the provided reaction if reaction in single compartment
                    - 'exchange' if reaction in multiple compartments

        (2) Case not a valid compartment:
                np.nan: 
                    'NaN' if one of the found compartments is not in :py:const:`~refinegems.utility.util.VALID_COMPARTMENTS`
    """
    
    metabs_from_reac = requests.get(BIGG_REACTIONS_URL + bigg_id, allow_redirects=False).json()['metabolites']
            
    comps = [comp_dict.get('compartment_bigg_id') for comp_dict in metabs_from_reac]  # Get all compartments for reaction
    contained_in_compartments = [(comp in VALID_COMPARTMENTS.keys()) for comp in comps]  # Get True for correct compartment        
    if not all(contained_in_compartments):  # At least one compartment not correct
        return np.nan
    else:  # All compartments correct
        if len(set(comps)) == 1:  # Set of found compartments of reaction = 1: Reaction happens in one compartment
            return comps[0]
        else:  # Not so important but do not remove reaction as reaction in correct compartments
            return 'exchange'  # Probably exchange reaction

def get_BiGG_metabs_annot_via_dbid(metabolite:cobra.Metabolite, 
                                        id:str, dbcol:str, 
                                        compartment:str = 'c') -> None:
    """Search for a BiGG ID and add it to a metabolite annotation. 
    The search is based on a column name of the BiGG metabolite table
    and an ID to search for. Additionally, using the given compartment name, 
    the found IDs are filtered for matching compartments.

    Args:
        - metabolite (cobra.Metabolite): 
            The metabolite. Needs to a a COBRApy Metabolte object.
        - id (str): 
            The ID to search for in the database.
        - dbcol (str): 
            Name of the column of the database to check the ID against.
        - compartment (str, optional): 
            The compartment name. Needs to be a valid BiGG compartment ID. 
            Defaults to 'c'.
    """
    # check, if a BiGG annotation is given
    if not 'bigg.metabolite' in metabolite.annotation.keys():
        # seach the database for the found id
        bigg_search = load_a_table_from_database(
            f'SELECT * FROM bigg_metabolites WHERE \'{dbcol}\' = \'{id}\'',
            query=True)
        if len(bigg_search) > 0:
            # check, if the matches also match the compartment
            metabolite.annotation['bigg.metabolite'] = [_.rsplit('_',1)[0] for _ in bigg_search['id'].tolist() if _.endswith(f'_{compartment}')]
            # add final matches to the annotations of the metabolite
            if len(metabolite.annotation['bigg.metabolite']) == 0:
                metabolite.annotation.pop('bigg.metabolite')
 
 
def add_annotations_from_BiGG_metabs(metabolite:cobra.Metabolite) -> None:
    """Check a cobra.metabolite for bigg.metabolite annotations. If they exists, 
    search for more annotations in the BiGG database and add them to the metabolite.

    Args:
        - metabolite (cobra.Metabolite): 
            The metabolite object.
    """
    if 'bigg.metabolite' in metabolite.annotation.keys():
        bigg_information = load_a_table_from_database(
            'SELECT * FROM bigg_metabolites WHERE id = \'' + f'\' OR id = \''.join(metabolite.annotation['bigg.metabolite']) + '\'',
            query=True)
        db_id_bigg = {'BioCyc':'biocyc', 'MetaNetX (MNX) Chemical':'metanetx.chemical','SEED Compound':'seed.compound','CHEBI':'chebi', 'KEGG Compound':'kegg.compound'}
        for db in db_id_bigg:
            info = list(set(bigg_information[db].dropna().to_list()))
            if len(info) > 0:
                info = ','.join(info)
                info = [x.strip() for x in info.split(',')] # make sure all entries are a separate list object
                if db_id_bigg[db] in metabolite.annotation.keys():
                    metabolite.annotation[db_id_bigg[db]] = list(set(info + metabolite.annotation[db_id_bigg[db]]))
                else:
                    metabolite.annotation[db_id_bigg[db]] = info


def _add_annotations_from_bigg_reac_row(row:pd.Series, reac:cobra.Reaction) -> None:
    """Given a row of the BiGG reaction database table and a cobra.Reaction object,
    extend the annotation of the latter with the information of the former.

    Args:
        - row (pd.Series): 
            The row of the database table.
        - reac (cobra.Reaction): 
            The reaction object.
    """
    
    dbnames = {'RHEA':'rhea','BioCyc':'biocyc','MetaNetX (MNX) Equation':'metanetx.reaction','EC Number':'ec-code'}
    for dbname,dbprefix in dbnames.items():
        if row[dbname]:
            ids_to_add = row[dbname].split(',')
            if dbprefix in reac.annotation.keys():
                reac.annotation[dbprefix] = list(set(reac.annotation[dbprefix]).union(set(ids_to_add)))
            else:
                reac.annotation[dbprefix] = ids_to_add


# @TEST
# @NOTE : A lot of warnings
def keep_only_bigg_reactions_in_certain_compartments(complete_df: pd.DataFrame) -> pd.DataFrame:
    """Extracts all possible BiGG ID variations from database for a BiGG reaction ID, gets the metabolite compartments 
    & returns table containing only reactions which happen in one of the provided compartments
        
    Args:
        - complete_df (pd.DataFrame): 
            Table containing at least the column 'bigg_id'.
        
    Returns:
        pd.DataFrame: 
            Table containing reactions & their compartments
    """
    tqdm.pandas()
    
    # (1) Find all occurrencs of a BiGG reaction ID in bigg_reactions table in database
    def get_all_similar_bigg_ids(bigg_id_in: str) -> list[str]:
        
        if '_' in bigg_id_in: bigg_id = re.split('_([a-zA-Z]|[0-9])$', bigg_id_in)[0]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): bigg_id = bigg_id_in[:-1]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): bigg_id = bigg_id_in[:-2]
        else: bigg_id = bigg_id_in
        
        query = f"SELECT id, INSTR(id, '{bigg_id}') bi FROM bigg_reactions WHERE bi > 0"
        result = con.execute(query).fetchall()
        result = [result_tuple[0] for result_tuple in result] if result else [bigg_id_in]
        result = [res for res in result if compare_bigg_ids(bigg_id, res)]
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
            for out in tqdm(pool.imap(get_reaction_compartment_from_bigg, complete_df.loc[:, "bigg_id"], chunksize=20), total=len(complete_df)):
                results.append(out)

        return results

    print(complete_df.columns)
    
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
    print(f'Getting all IDs with correct compartment {VALID_COMPARTMENTS.keys()}...')
    results = multi_get_reaction_compartment(complete_df)
    complete_df["compartment"] = results

    # complete_df.progress_apply(get_reaction_compartment, axis=1)  # (2)

    # (3) Remove reactions with compartment = NaN
    complete_df.dropna(subset=['compartment'], inplace=True)

    return complete_df




# @DISCUSSION
# -----------
def compare_bigg_ids(id1: str, id2: str) -> bool:
    """Compares two BiGG strings/IDs & Returns True if one BiGG ID matches most of the other

    Args:
        - id1 (str): 
            ID 1
        - id2 (str): 
            ID 2

    Returns:
        bool: 
            Indicates if most of one string contained in the other
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


# @TEST
def get_bigg_db_mapping(map_to:str='BioCyc', metabolites:bool=True) -> pd.DataFrame:
    """Download a mapping of BiGG IDs to a specified database.

    Args:
        map_to (str, optional): 
            Name of the database to map to. 
            Ideally a column of the table in the database, 
            but SEED, KEGG and BioCyc are valid as well. 
            Defaults to 'BioCyc'.
        metabolites (bool, optional): 
            Flag to map reaction (False) or metabolite (True) IDs. 
            Defaults to True.

    Raises:
        - KeyError: Given database name not found in database. Cannot perform mapping.

    Returns:
        pd.DataFrame: 
            The mapping as a table.
    """

    # adjust name to map to if necessary
    reac_or_comp = 'Compound' if metabolites else 'Reaction'
    table_name = 'bigg_metabolites' if metabolites else 'bigg_reactions'
    if map_to in ['SEED','KEGG','Reactome']:
        map_to =  ' '.join([map_to, reac_or_comp])
    
    # download BiGG tables from database
    # ----------------------------------
    # build connection to DB
    connection = sqlite3.connect(PATH_TO_DB)
    cursor = connection.cursor()

    # retrieve only mappings to a specific database
    result = cursor.execute('SELECT 1 FROM PRAGMA_TABLE_INFO(?) WHERE name = ?',(table_name,map_to))
    possible_db = result.fetchone()
    if possible_db:
        query = f'SELECT * FROM {table_name} WHERE {map_to} IS NOT NULL'
    else:
        raise KeyError('Given database name not found in database. Cannot perform mapping.')

    # actually load data
    data = load_a_table_from_database(query)
    data = data.explode(map_to, ignore_index=True)

    # reduce columns to mapping only
    data = data[['id',map_to]]
    data.rename(columns={'id':'bigg_id'}, inplace=True)

    # filter for compartment in case of reactions
    if not metabolites:
        data = keep_only_bigg_reactions_in_certain_compartments(data)

    # close connection to database
    connection.close()
    
    return data


# ChEBI
# -----
"""Primarily used for curating the media database. Obtain information about 
metabolites (table format)."""

# @DISCUSSION: Add to metabolite handling for `.entities.build_metabolite_bigg`?
def add_charges_chemical_formulae_to_metabs(missing_metabs: pd.DataFrame) -> pd.DataFrame:
   """Adds charges & chemical formulae from CHEBI/BiGG to the provided dataframe

   Args:
      - missing_metabs (pd.DataFrame): 
         Table containing metabolites & the respective CHEBI & BiGG IDs
         
   Returns:
      pd.DataFrame: 
         Input table extended with the charges & chemical formulas obtained from CHEBI/BiGG
   """
   
   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id = str(row.get('bigg_id'))
      charge = None
      if chebi_id:  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         return chebi_entity.get_charge()
      elif bigg_id != 'nan':  # Get charge from BiGG if no ChEBI ID available
         try:
            charge = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['charges'][0]  # Take first charge
         except ValueError:
            pass   
         # If no charge was found, charge=0
         return charge if charge else 0
   
   # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
   def find_formula(row: pd.Series) -> str:
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id, chem_form = str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id: # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except (ValueError, IndexError) as e:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else '*'
      return chem_formula
   
   missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
   missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
   missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   
   return missing_metabs

def add_info_from_ChEBI_BiGG(missing_metabs: pd.DataFrame, charge=True, formula=True, iupac=True) -> pd.DataFrame:
   """Adds information from CHEBI/BiGG to the provided dataframe.

   The following informations can be added:

   - charge
   - formula
   - iupac (name)

   Args:
      - missing_metabs (pd.DataFrame): 
         Table containing metabolites & the respective ChEBI & BiGG IDs
         
   Returns:
      pd.DataFrame: 
         Input table extended with the charges & chemical formulas obtained from ChEBI/BiGG.
   """

   # check if a row contains a ChEBI ID, take the first and make sure its in the format: CHEBI:234567
   def get_chebi_id(row: pd.Series) -> str:

      chebi = row.get('ChEBI')
      if pd.isnull(chebi):
         return None
      elif type(chebi) == str: # check for mulitple entries (assuming first one is the best fitting one)
         chebi = chebi.split(',')[0]
         if 'CHEBI' in chebi:
            return chebi
         else:
            return 'CHEBI:' + chebi
      else:
         return 'CHEBI:' + str(chebi)

   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id = get_chebi_id(row)
      bigg_id = str(row.get('bigg_id'))
      charge = None
      if chebi_id:  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         return chebi_entity.get_charge()
      elif bigg_id != 'nan':  # Get charge from BiGG if no ChEBI ID available
         try:
            charge = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['charges'][0]  # Take first charge
         except ValueError:
            pass   
         # If no charge was found, charge=0
         return charge if charge else 0
    
    # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
   def find_formula(row: pd.Series) -> str:
      chebi_id = get_chebi_id(row)
      bigg_id, chem_form = str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id: # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except ValueError:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else 'No formula'
      return chem_formula
   
   # using the ChEBI ID, retrieve the IUPAC name from the ChEBI database.
   def find_iupac(row: pd.Series) -> str:

      chebi_id = get_chebi_id(row)
      if chebi_id:
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         # only take the IUPAC names
         iupac_names = []
         for name in chebi_entity.get_names():
            print(name)
            if name.get_source() == 'IUPAC':
               iupac_names.append(name.get_name())
               break

         if len(iupac_names) > 0:
            iupac_name = ', '.join(iupac_names)
         else:
            iupac_name = None
      else:
         iupac_name = None

      return iupac_name

   if charge:
      missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   if formula:
      missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
      missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
      missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   if iupac:
      missing_metabs['ChEBI_IUPAC'] = missing_metabs.apply(find_iupac, axis=1)
   
   return missing_metabs


# KEGG
# ----
"""Retrieve and parse information from the KEGG database. 
"""

def get_kegg_genes(organismid: str) -> pd.DataFrame:
    """Extracts list of genes from KEGG given an organism

    Args:
        - organismid (str): 
            KEGG ID of organism which the model is based on

    Returns:
        pd.DataFrame: 
            Table of all genes denoted in KEGG for the organism
    """
    k = KEGG()
    gene_list = k.list(organismid)

    return pd.read_table(io.StringIO(gene_list), header=None)


def parse_KEGG_gene(locus_tag:str) -> dict:
    """Based on a locus tag, fetch the corresponding KEGG entry and
    parse it into a dictionary containing the following information (if available):
    
    - ec-code
    - orthology
    - references

    Args:
        - locus_tag (str): 
            The locus in the format <organism_id>:<locus_tag>

    Returns:
        dict: 
            The collected information.
    """
    
    gene_info = dict()
    gene_info['orgid:locus'] = locus_tag
    
    # retireve KEGG gene entry 
    try: 
        gene_entry = list(Gene.parse(REST.kegg_get(locus_tag)))[0]
    except Exception as e:
        # @TODO : warning / logging
        gene_entry = None
    
    # skip, if no entry found
    if not gene_entry:
        gene_info['ec-code'] = None
        return gene_info
    
    # extract orthology and ec-code
    if len(gene_entry.orthology) > 0:
        # get KEGG orthology ID
        kegg_orthology = [_[0] for _ in gene_entry.orthology]
        gene_info['kegg.orthology'] = kegg_orthology
        # get EC number
        ec_numbers = [re.search('(?<=EC:).*(?=\])',orth[1]).group(0) for orth in gene_entry.orthology if re.search('(?<=EC:).*(?=\])',orth[1])]
        if isinstance(ec_numbers,list) and len(ec_numbers) > 0:
            gene_info['ec-code'] = [ec for ec_str in ec_numbers for ec in ec_str.split(' ')]
            
    if not 'ec-code' in gene_info.keys():
        gene_info['ec-code'] = None
        
    # get more information about connections to other databases
    if len(gene_entry.dblinks) > 0:
        for dbname,ids in gene_entry.dblinks:
            conform_dbname = re.sub(pattern='(NCBI)(.*)(ID$)', repl='\\1\\2',string=dbname) # Remove ID if NCBI in name
            conform_dbname = re.sub('[^\w]','',conform_dbname) # remove special signs except underscore
            conform_dbname = conform_dbname.lower() # make lower case
            gene_info[conform_dbname] = ids
            
    return gene_info
            
    
def parse_KEGG_ec(ec:str) -> dict:
    """Based on an EC number, fetch the corresponding KEGG entry and
    parse it into a dictionary containing the following information (if available):
    
    - ec-code
    - id (kegg.reference)
    - equation
    - reference 
    - pathway

    Args:
        - ec (str): 
            The EC number in the format 'x.x.x.x'

    Returns:
        dict: 
            The collected information about the KEGG entry.
    """
    
    ec_info = dict()
    ec_info['ec-code'] = ec
    
    # retrieve KEGG entry
    # @TODO : add time restraint and tries for time out only
    try:
        ec_entry = list(Enzyme.parse(REST.kegg_get(ec)))[0]
    except Exception as e:
        # @TODO logging / warning
        ec_entry = None
        ec_info['id'] = None
        ec_info['equation'] = None
        ec_info['reference'] = None
        return ec_info
    
    # retrieve reaction information from entry 
    rn_numbers = [re.search('(?<=RN:).*(?=\])',reac).group(0) for reac in ec_entry.reaction if re.search('(?<=RN:).*(?=\])',reac)]
    if len(rn_numbers) > 0:
        ec_info['id'] = [_.split(' ') for _ in rn_numbers]
        ec_info['equation'] = None
    else:
        ec_info['id'] = None
        ec_info['equation'] = ec_entry.reaction
        
    # retrieve database links from entry
    refs = dict()
    # orthology not possible with biopython
    if len(ec_entry.pathway) > 0:
        refs['kegg.pathway'] = [_[1] for _ in ec_entry.pathway]
    # @TODO extend as needed
    if len(ec_entry.dblinks) > 0:
        for dbname, ids in ec_entry.dblinks:
            if 'BRENDA' in dbname:
                refs['brenda'] = ids
            if 'CAS' == dbname:
                refs['cas'] = ids
    ec_info['reference'] = refs
    
    return ec_info


# @TODO open issues
def kegg_reaction_parser(rn_id:str) -> dict: 
    """Get the entry of a KEGG reaction ID and 
    parse the information into a dictionary.

    Args:
        - rn_id (str): 
            A reaction ID existing in KEGG.

    Returns:
        dict: 
            The KEGG entry information as a dictionary.
    """

    # get KEGG reaction entry
    try:
        kegg_reac = REST.kegg_get(F'rn:{rn_id}')
        kegg_reac = kegg_reac.read()
    except Exception as e:
        # @TODO
        return None

    # parse the entry for necessary information
    features = {}
    db_entries = []
    pathways = []
    rc = []
    references = {'kegg.reaction':rn_id}
    collect = False
    for line in kegg_reac.split('\n'):
        if line:
            if line.startswith('NAME'):
                features['name'] = line.replace('NAME','',1).strip()
            elif line.startswith('EQUATION'):
                features['equation'] = line.replace('EQUATION','',1).strip()
            elif line.startswith('ENZYME'):
                references['ec-code'] = line.replace('ENZYME','',1).strip()
            elif line.startswith('RCLASS'):
                rc.append(line.replace('RCLASS','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('PATHWAY'):
                pathways.append(line.replace('PATHWAY','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('DBLINKS'):
                db_entries.append(line.replace('DBLINKS','',1).strip())
                collect = True
            elif collect == True and line[0] != '/':
                if len(db_entries) == 0:
                    if line[0].isupper():
                        collect = False
                    else:
                        line = line.strip()
                        if line.startswith('RC'):
                            rc.append(line.split(' ')[0])
                        else:
                            pathways.append(line.split(' ')[0])
                else:
                    db_entries.append(line.strip())
            else:
                continue

    # parse references
    for entry in db_entries:
        db, identifier = entry.split(':')
        db = db.strip().lower()
        if db in references:
            references[db] = references[db].append(identifier)
        else:
            references[db] = [identifier.strip()]
    if len(pathways) > 0:
        references['kegg.pathway'] = pathways
    if len(rc) > 0:
        references['kegg.rclass'] = rc
    if len(references) > 0:
        features['db'] = references

    return features


# ModelSEED
# ---------
""" Reports mismatches in charges and formulae based on ModelSEED

Extracts ModelSEED data from a given tsv file, extracts all metabolites from a given model. Both lists of metabolites are compared by charge and formula.
"""

def get_modelseed_compounds() -> pd.DataFrame:
    """Extracts compounds from ModelSEED which have BiGG Ids

    Returns:
        pd.DataFrame: 
            Table containing ModelSEED data
    """
    # Get only rows where BiGG is contained
    com = load_a_table_from_database("SELECT *, INSTR(aliases, 'BiGG:') bigg FROM modelseed_compounds WHERE bigg > 0")

    def get_bigg_ids(aliases):
        try:
            aliases_list = aliases.split('|')
            bigg = [x[6:] for x in aliases_list if re.search('BiGG: .*', x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com['BiGG'] = com.apply(lambda row: get_bigg_ids(row['aliases']), axis=1)

    return com.loc[:, ['id', 'name', 'formula', 'mass', 'charge', 'BiGG']]


def get_model_charges(model: cobraModel) -> pd.DataFrame:
    """Extracts all metabolites from model

    Args:
        - model (cobraModel): 
            Model loaded with COBRApy

    Returns:
        pd.DataFrame: 
            Table containing charges and formulae of model metabolites
    """
    charges = {}
    for metab in model.metabolites:
        charges[metab.id[:-2]] = [metab.charge, metab.formula]

    df_charges = pd.DataFrame.from_dict(
        charges,
        orient='index',
        columns=[
            'charge_model',
            'formula_model']).reset_index().rename(
        columns={
            'index': 'BiGG'})

    return df_charges


def get_modelseed_charges(modelseed_compounds: pd.DataFrame) -> pd.DataFrame:
    """Extract table with BiGG, charges and formulae

    Args:
        - modelseed_compounds (pd.DataFrame): 
            ModelSEED data. Output from :py:func:`~refinegems.utility.db_access.get_modelseed_compounds`.

    Returns:
        pd.DataFrame: 
            Table containing charges and formulae of ModelSEED metabolites
    """
    modelseed_compounds = modelseed_compounds.loc[:, ['charge', 'BiGG', 'formula']].rename(
        columns={'charge': 'charge_modelseed', 'formula': 'formula_modelseed'})
    lst_col = 'BiGG'
    x = modelseed_compounds.assign(
        **{lst_col: modelseed_compounds[lst_col].str.split(';')})
    df_ms = pd.DataFrame({col: np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference(
        [lst_col])}).assign(**{lst_col: np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return df_ms


def compare_model_modelseed(model_charges: pd.DataFrame, modelseed_charges: pd.DataFrame) -> pd.DataFrame:
    """Compares tables with charges / formulae from model & modelseed

    Args:
        - model_charges (pd.DataFrame): 
            Charges and formulae of model metabolites. Output of :py:func:`~refinegems.utility.db_access.get_model_charges`.
        - modelseed_charges (pd.DataFrame): 
            Charges and formulae of ModelSEED metabolites. Output of :py:func:`~refinegems.utility.db_access.get_modelseed_charges`.

    Returns:
        pd.DataFrame: 
            Table containing info whether charges / formulae match
    """
    df_comp = pd.merge(model_charges, modelseed_charges, on='BiGG', how='left')

    def f(x):
        return True if float(
            x['charge_model']) == x['charge_modelseed'] else False

    def g(x):
        return True if x['formula_model'] == x['formula_modelseed'] else False

    df_comp['charge_match'] = df_comp.apply(f, axis=1)
    df_comp['formula_match'] = df_comp.apply(g, axis=1)

    return df_comp


def get_charge_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with charge mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame): 
            Charge and formula mismatches. Output from :py:func:`~refinegems.utility.db_access.compare_model_modelseed`.

    Returns:
        pd.DataFrame: 
            Table containing metabolites with charge mismatch
    """
    return df_comp.loc[~df_comp['charge_match']].dropna(
        subset=['charge_modelseed'])


def get_formula_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with formula mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame): 
            Charge and formula mismatches. Output from :py:func:`~refinegems.utility.db_access.compare_model_modelseed`.

    Returns:
        pd.DataFrame: 
            Table containing metabolites with formula mismatch
    """
    return df_comp.loc[~df_comp['formula_match']].dropna(
        subset=['formula_modelseed'])


def get_compared_formulae(formula_mismatch: pd.DataFrame) -> pd.DataFrame:
    """Compare formula by atom pattern

    Args:
        formula_mismatch (pd.DataFrame): 
            Table with column containing atom comparison. Output from :py:func:`~refinegems.utility.db_access.get_formula_mismatch`.

    Returns:
        pd.DataFrame: 
            table containing metabolites with formula mismatch
    """

    def formula_comparison(f1, f2):
        # from Jan Leusch
        formula_pattern = "[A-Z][a-z]?\\d*"
        atom_pattern = '[A-Z][a-z]?'
        atom_number_pattern = '\\d+'
        difference = {}
        f1_dict = {}
        f2_dict = {}
        f1 = re.findall(formula_pattern, f1)
        f2 = re.findall(formula_pattern, f2)
        for p in f1:
            key = re.findall(atom_pattern, p)[0]
            value = re.findall(atom_number_pattern, p)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f1_dict[key] = value
        for q in f2:
            key = re.findall(atom_pattern, q)[0]
            value = re.findall(atom_number_pattern, q)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f2_dict[key] = value
        difference = f1_dict
        if not f2_dict:
            return difference
        else:
            for q in f2_dict:
                if q in difference:
                    difference[q] -= f2_dict[q]
                else:
                    difference[q] = -f2_dict[q]
        for a in list(difference):
            if difference[a] == 0:
                difference.pop(a)
        return difference

    formula_mismatch['formula_comparison'] = formula_mismatch.apply(
        lambda row: formula_comparison(
            row['formula_model'], row['formula_modelseed']), axis=1)

    return formula_mismatch


def compare_to_modelseed(model: cobraModel) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Executes all steps to compare model metabolites to ModelSEED metabolites

    Args:
        - model (cobraModel): 
            Model loaded with COBRApy

    Returns:
        tuple: 
            Tables with charge (1) & formula (2) mismatches

            (1) pd.DataFrame: Table with charge mismatches 
            (2) pd.DataFrame: Table with formula mismatches
    """
    ms_comp = get_modelseed_compounds()
    model_charges = get_model_charges(model)
    modelseed_charges = get_modelseed_charges(ms_comp)
    df_comp = compare_model_modelseed(model_charges, modelseed_charges)
    charge_mismatch = get_charge_mismatch(df_comp)
    formula_mismatch = get_formula_mismatch(df_comp)
    formula_comp = get_compared_formulae(formula_mismatch)
    return charge_mismatch, formula_comp



# NCBI
# ----
"""Connect to NBCI (using Biopython/Entrez) to extract information
like EC number of an NBCI protein accession number or information about locus tags.
"""

def search_ncbi_for_gpr(locus: str) -> str:
    """Fetches protein name from NCBI

    Args:
        - locus (str): 
            NCBI compatible locus_tag

    Returns:
        str: 
            Protein name|description
    """
    handle = Entrez.efetch(
        db="protein",
        id=locus,
        rettype="gbwithparts",
        retmode='text')
    records = SeqIO.parse(handle, "gb")

    for i, record in enumerate(records):
        if (locus[0] == 'W'):
            return record.description, locus
        else:
            for feature in record.features:
                if feature.type == "CDS":
                    return record.description, feature.qualifiers["locus_tag"][0]


# fetching the EC number (if possible) from NCBI 
# based on an NCBI protein ID (accession version number)
# @TODO logging
def get_ec_from_ncbi(mail:str,ncbiprot:str) -> str|None:
    """Based on a NCBI protein accession number, try and fetch the 
    EC number from NCBI.

    Args:
        - mail (str): 
            User's mail address for the NCBI ENtrez tool.
        - ncbiprot (str): 
            The NCBI protein accession number.

    Returns:
        (1) Case: fetching successful 
                str: 
                    The EC number associated with the protein ID based on NCBI.
        
        (2) Case: fetching unsuccessful
                None:
                    Nothing to return
     """
    try: 
        Entrez.email = mail
        handle = Entrez.efetch(db='protein', id=ncbiprot, rettype='gpc', retmode='xml')
        for feature_dict in xmltodict.parse(handle)['INSDSet']['INSDSeq']['INSDSeq_feature-table']['INSDFeature']:
            if isinstance(feature_dict['INSDFeature_quals']['INSDQualifier'], list):
                for qual in feature_dict['INSDFeature_quals']['INSDQualifier']:
                    if qual['INSDQualifier_name'] == 'EC_number':
                        return qual['INSDQualifier_value']
    except Exception as e:
        # @TODO : logging / warning etc. ?
        return None
  
  
# Uniprot
# -------
"""Currently mainly used for gapfilling, map your proteins against EC numbers 
using the SwissProt database"""

def map_dmnd_res_to_sp_ec_brenda(dmnd_results:pd.DataFrame, 
                                 swissprot_mapping_path:str) -> pd.DataFrame:
    """Map the results of a DIAMOND BLASTp run (filtered, see 
    :py:func:`~refinegems.curation.db_access.db.filter_DIAMOND_blastp_results`)

    Args:
        - dmnd_results (pd.DataFrame): 
            The results of the DIAMOND run.
        - swissprot_mapping_path (str): 
            The path to the SwissProt mapping file (IDs against BRENDA and EC,
            for information on how to get them, refer to :py:func:`~refinegems.utility.set_up.download_url`)

    Returns:
        pd.DataFrame: 
            The resulting mapping (no duplicates).
    """
    
    def _combine_EC_BRENDA(brenda:str,ec:str) -> list:
        """Helper function for combining the entries for BRENDA and EC number 
        for one ID (SwissProt mapping file) into a list of EC numbers.

        Args:
            - brenda (str): 
                The BRENDA column value.
            - ec (str): 
                The EC number column value.

        Returns:
            list: 
                The list with the EC numbers.
        """
    
        nums_list = []
        # get brenda ids
        if brenda and not pd.isna(brenda):
            for num in brenda.split(';'):
                nums_list.append(num.strip())

        # get ec numbers
        if ec and not pd.isna(ec):
            for num in ec.split(';'):
                nums_list.append(num.strip())
            
        # remove duplicates
        nums_list = list(set(nums_list))
        # filter non-complete EC numbers : X.X.X.X => 7 or more characters
        nums_list = [_ for _ in nums_list if len(_) >= 7]
            
        return nums_list
    
    # load the SwissProt mapping file
    swissprot_mapping = pd.read_csv(swissprot_mapping_path, sep='\t')
    swissprot_mapping.dropna(subset=['BRENDA','EC number'], how='all',inplace=True)
    
    # extract SwissProts IDs from subject_ID
    dmnd_results.columns = ['locus_tag','UniProt']
    dmnd_results['UniProt'] = dmnd_results['UniProt'].apply(lambda x: x.split('|')[1])
    
    # match
    dmnd_results = dmnd_results.merge(swissprot_mapping, left_on='UniProt', right_on='Entry', how='left')
    dmnd_results.drop('Entry', axis=1, inplace=True)
    dmnd_results['ec-code'] = dmnd_results.apply(lambda x: _combine_EC_BRENDA(x['BRENDA'],x['EC number']), axis=1)
    dmnd_results.drop(['BRENDA','EC number'], axis=1, inplace=True)
    dmnd_results = dmnd_results.explode('ec-code')
    dmnd_results.drop_duplicates(inplace=True)
    
    return dmnd_results

def get_ec_via_swissprot(fasta:str, db:str, missing_genes:pd.DataFrame, 
                         swissprot_mapping_file:str,
                         outdir:str=None,
                         sens:Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive']='more-sensitive',
                         cov:float=95.0,
                         t:int=2, pid:float=90.0) -> pd.DataFrame:
    """Based on a protein FASTA and a missing genes tables, mapped them to EC numbers 
    using a Swissprot DIAMOND database and a SwissProt mapping file (see :py:func:`~refinegems.utility.set_up.download_url`
    on how to download the needed files).

    Args:
        - fasta (str): 
            Path to the FASTA protein file.
        - db (str): 
            Path to the DIAMOND database (SwissProt).
        - missing_genes (pd.DataFrame): 
            The table of missing genes.
        - swissprot_mapping_file (str): 
            Path to the SwissProt mapping file.
        - outdir (str, optional): 
            Path to a directory to write the output to. 
            Defaults to None.
        - sens (Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], optional): 
            Sensitivity mode of DIAMOND blastp. Defaults to 'more-sensitive'.
        - cov (float, optional): 
            Coverage threshold for DIAMOND blastp. Defaults to 95.0.
        - t (int, optional): 
            Number of threads to use for DIAMOND blastp. Defaults to 2.
        - pid (float, optional): 
            Percentage identity value to use as a cutoff for the results
            of the DIAMOND blastp run. Defaults to 90.0.

    Returns:
        pd.DataFrame: 
            The missing genes table extended by the mapping to an EC number, 
            if successful.
    """
    
    # Step 1: Make a FASTA out of the missing genes
    miss_fasta = create_missing_genes_protein_fasta(fasta,outdir,missing_genes)
    # Step 2: Run DIAMOND
    #         blastp mode against SwissProt DB
    blast_path = run_DIAMOND_blastp(miss_fasta, db, 
                                    sensitivity=sens,
                                    coverage=cov,
                                    threads=t,
                                    outdir=outdir)
    # Step 3: filter DIAMOND hits
    dmnd_res = filter_DIAMOND_blastp_results(blast_path, pid)
    # Step 4: map to Swissprot mapping file
    mapped_res = map_dmnd_res_to_sp_ec_brenda(dmnd_res, swissprot_mapping_file)
    # Step 5: Aggregate UniProt IDs for unique combinations of 
    #         EC numbers and locus tags
    mapped_res = mapped_res.groupby(['locus_tag','ec-code']).agg({'UniProt': lambda x: x.tolist()}).reset_index()

    return mapped_res

