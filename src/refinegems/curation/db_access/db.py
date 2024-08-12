#!/usr/bin/env python
"""This module contains functions usable with multiple databases or those drawing
information from multiple databases and therefore cannot be separated into the other 
submodules. Additionally, this module contains function for databases not covered by
the other modules.
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune and Tobias Fehrenbach"

################################################################################
# requirements
################################################################################

import cobra
import numpy as np
import pandas as pd
import re
import requests
import sqlite3
import subprocess
import xmltodict

from Bio import Entrez
from multiprocessing import Pool
from pathlib import Path
from ratelimit import limits, sleep_and_retry
from tqdm import tqdm
from typing import Literal

tqdm.pandas()
pd.options.mode.chained_assignment = None # suppresses the pandas SettingWithCopyWarning; comment out before developing!!

from ...utility.databases import PATH_TO_DB
from ...utility.entities import VALID_COMPARTMENTS
from ...utility.io import load_a_table_from_database, create_missing_genes_protein_fasta

################################################################################
# variables
################################################################################

ALL_BIGG_COMPARTMENTS_ONE_LETTER = ('c', 'e', 'p', 'm', 'x', 'r', 'v', 'n', 'g', 'u', 'l', 'h', 'f', 's', 'i', 'w', 'y') #: :meta: 
ALL_BIGG_COMPARTMENTS_TWO_LETTER = ('im', 'cx', 'um', 'cm', 'mm') #: :meta: 
BIGG_REACTIONS_URL = 'http://bigg.ucsd.edu/api/v2/universal/reactions/' #: :meta: 
BIGG_METABOLITES_URL = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/' #: :meta: 

################################################################################
# functions
################################################################################
      
def compare_ids(id1: str, id2: str) -> bool:
    """Compares two strings/IDs & Returns True if one string matches most of the other

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
 

# BiGG
# ----

@sleep_and_retry
@limits(calls=10, period=1)
def get_reaction_compartment(bigg_id: str) -> str:
    """Retrieves the compatment(s) a reaction is hapening in from the BiGG reaction identifier
        via the metabolites

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
                    'NaN' if one of the found compartments is not in VALID_COMPARTMENTS.keys()
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

# @TEST
# @NOTE : A lot of warnings
def keep_only_reactions_in_certain_compartments(complete_df: pd.DataFrame) -> pd.DataFrame:
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
        data = keep_only_reactions_in_certain_compartments(data)

    # close connection to database
    connection.close()
    
    return data

 
# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_bigg_model(complete_df: pd.DataFrame, model_entities: pd.DataFrame, metabolites: bool=False) -> pd.DataFrame:
    """Compares missing entities obtained through genes extracted via KEGG/BioCyc to entities in the model
        Needed to back check previous comparisons.

    Args:
        - complete_df (pd.DataFrame): 
            Table that contains KEGG/BioCyc Id, BiGG Id & more
        - model_entities (pd.DataFrame): 
            BiGG Ids of entities in the model 
        - metabolites (bool): 
            True if names of metabolites should be added, otherwise false

    Returns:
        pd.DataFrame: 
            Table containing entities present in KEGG/BioCyc but not in the model
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
            return compartment if compartment in VALID_COMPARTMENTS.keys() else np.nan  # To filter the incorrect compartments out
        
        entities_missing_in_model['compartment'] = entities_missing_in_model.apply(
            lambda row: get_compartment_from_id(row['bigg_id']), axis=1)
        entities_missing_in_model.dropna(subset=['compartment'], inplace=True)  # Drop all BiGG metabolite IDs which have no valid compartment
    
    return entities_missing_in_model


def add_stoichiometric_values_to_reacs(missing_reacs: pd.DataFrame) -> pd.DataFrame:
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
 
 
 # @TODO : name is an issue

 
# NCBI parser
# -----------

# fetching the EC number (if possible) from NCBI 
# based on an NCBI protein ID (accession version number)
def get_ec_from_ncbi(mail:str,ncbiprot:str):
    try: 
        Entrez.email = mail
        handle = Entrez.efetch(db='protein', id=ncbiprot, rettype='gpc', retmode='xml')
        for feature_dict in xmltodict.parse(handle)['INSDSet']['INSDSeq']['INSDSeq_feature-table']['INSDFeature']:
            if isinstance(feature_dict['INSDFeature_quals']['INSDQualifier'], list):
                for qual in feature_dict['INSDFeature_quals']['INSDQualifier']:
                    if qual['INSDQualifier_name'] == 'EC_number':
                        return qual['INSDQualifier_value']
    except Exception as e:
        # @TODO : logging / warning etc
        return None
    
    
# DIAMOND
# -------
# @ISSUE / @NOTE / @TODO
#   putting these in connections leads to an import error

def run_DIAMOND_blastp(fasta:str, db:str, 
                       sensitivity:Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive']='more-sensitive',
                       coverage:float=95.0,
                       threads:int=2,
                       outdir:str=None, outname:str='DIAMOND_blastp_res.tsv'):
    
    if outdir:
        outname = Path(outdir,'DIAMOND_blastp_res.tsv')
        logfile = Path(outdir,'log_DIAMOND_blastp.txt')
    else:
        outname = Path(outname)
        logfile = Path('log_DIAMOND_blastp.txt')
      
    # @TODO: test, if it works with different paths and their problems  
    # @TODO: write additional output to a logfile, not stderr
    subprocess.run([F'diamond blastp -d {db} -q {fasta} --{sensitivity} --query-cover {coverage} -p {int(threads)} -o {outname} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore 2> {logfile}'], shell=True)

    return outname

def filter_DIAMOND_blastp_results(blasttsv:str, pid_theshold:float=90.0):
    
    if pid_theshold > 100.0 or pid_theshold < 0.0:
        raise ValueError('PID threshold has to be between 0.0 and 100.0')
    
    # load diamond results
    diamond_results = pd.read_csv(blasttsv, sep='\t', header=None)
    diamond_results.columns = ['query_ID', 'subject_ID', 'PID', 'align_len', 'no_mismatch', 'no_gapopen', 'query_start', 'query_end', 'subject_start', 'subject_end','E-value','bitscore']
    # filter by PID
    diamond_results = diamond_results[diamond_results['PID']>=pid_theshold]
    # trim cols
    diamond_results = diamond_results[['query_ID','subject_ID']]
    
    return diamond_results

# Uniprot
# -------

def map_dmnd_res_to_sp_ec_brenda(dmnd_results:pd.DataFrame, 
                                 swissprot_mapping_path:str):
    
    def _combine_EC_BRENDA(brenda:str,ec:str):
    
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


