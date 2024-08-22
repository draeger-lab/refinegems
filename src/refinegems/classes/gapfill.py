#!/usr/bin/env python
# @TODO: Correct/Rewrite doc string
# @TEST - new runtimes 
"""Add more reactions, genes and more to a model based on different gap-filling methods.
    The `gapfill` module can be used either with KEGG were you only need the KEGG organism code or with BioCyc or with both (Options: 'KEGG', 'BioCyc', 'KEGG+BioCyc').
    For how to obtain the BioCyc tables look into the documentation under: 'Main modules' > 'Gap filling' > 'Automated gap filling'.
    
    Run times: 
    
        * 'KEGG': ~ 2h
        * 'BioCyc': ~ 45mins - 1h
        * 'KEGG+BioCyc': ~ 3 - 4h  
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune and Dr. Reihaneh Mostolizadeh"

############################################################################
# requirements
############################################################################

from abc import ABC, abstractmethod
from operator import index

import cobra
import io
import numpy as np
import pandas as pd
import re
import sqlite3
import warnings

from bioservices.kegg import KEGG
from itertools import chain
from libsbml import Model as libModel
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal, Union

from tqdm import tqdm
tqdm.pandas()

from ..utility.databases import PATH_TO_DB
from ..utility.db_access import get_ec_from_ncbi, get_ec_via_swissprot, parse_KEGG_gene, parse_KEGG_ec, add_stoichiometric_values_to_reacs_from_bigg
from ..utility.io import load_a_table_from_database, parse_gff_for_cds, load_model, write_model_to_file
from ..utility.entities import create_gp, create_gpr, build_reaction_bigg, build_reaction_kegg, build_reaction_mnx, isreaction_complete, extract_metabolites_from_reactions
from ..utility.util import VALID_COMPARTMENTS
from ..developement.decorators import *

##### Coloured text:
from colorama import init as colorama_init
from colorama import Fore

def coloured_example_text():
    colorama_init(autoreset=True)
    print(f'{Fore.RED}To use the KEGG comparison the specification of the organismid (KEGG organism code) is obligatory.\n' +
      'If there is no organism code available for your organism in KEGG but an entry for your organism exists in BioCyc, use the option \'BioCyc\'.\n' +
      'If no entry for your organism exists in KEGG and/or BioCyc, the gap analysis cannot be done.')


# @Note:
#   some reactions have @DEBUGGING 
#   enable these parts to shorten runtime during debugging (used to work on a subset, 
#   not on the whole input)

############################################################################
# variables
############################################################################


############################################################################
# functions
############################################################################

# @DISCUSSION: Keep for end user? Use somewhere else? Or deprecate?
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


# Mapping for BioCyc Reactions
# ----------------------------
def map_biocyc_to_reac(biocyc_reacs: pd.DataFrame) -> tuple[pd.DataFrame, dict[str: int]]:

    statistics = {
        'mapped2MNX': 0, 
        'mapped2BiGG': 0, 
        'remaining_unmapped': 0
    }

    # Drop NaNs in relevant columns
    biocyc_reacs.dropna(subset=['id', 'ec-code'], inplace=True)

    def _map_to_mnx(biocyc_reacs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Maps Biocyc IDs in a table to the corresponding MetaNetX IDs & 
        Adds information obtained via the MetaNetX IDs to the resulting table

        Args:
            - biocyc_reacs (pd.DataFrame): 
                Table containing BioCyc IDs

        Returns:
            pd.DataFrame: 
                Table with BioCyc IDs mapped to MetaNetX IDs & more information 
                obtained from MetaNetX
        """

        # Step 1: Mapping & Obtain information
        # ------------------------------------
        # Mapping to MetaNetX + addition of more information on MetaNetX 
        # Reactions
        # Get MetaNetX BioCyc
        mnx2biocyc_reacs = load_a_table_from_database(
            '''
            SELECT x.source, x.id, p.mnx_equation, p.reference, p.\'ec-code\' 
            FROM mnx_reac_xref x INNER JOIN mnx_reac_prop p 
            USING(id) 
            WHERE x.source LIKE \'%metacyc%\' 
            OR x.source LIKE \'%biocyc%\'
            '''
            )

        # Change column names to fit input table
        mnx2biocyc_reacs.rename({'id': 'MetaNetX', 'source': 'id'}, axis=1, inplace=True)

        # Transform id column to only contain BioCyc/MetaCyc IDs
        mnx2biocyc_reacs['id'] = mnx2biocyc_reacs['id'].str.split(':', 1).str[-1]

        # Crop table to contain BiGG IDs set per BioCyc ID
        mnx_as_list = mnx2biocyc_reacs.groupby('id')['MetaNetX'].apply(set).reset_index(name='MetaNetX')
        mnx2biocyc_reacs.drop('MetaNetX', axis=1, inplace=True)
        mnx2biocyc_reacs = mnx_as_list.merge(mnx2biocyc_reacs, on='id')

        # Drop duplicates to get the unique BioCyc IDs
        mnx2biocyc_reacs.drop_duplicates(subset='id', inplace=True)

        print(mnx2biocyc_reacs[mnx2biocyc_reacs['id'].duplicated()])

        # Merge mnx2biocyc_reacs with biocyc_reacs to get new information
        biocyc_reacs = mnx2biocyc_reacs.merge(biocyc_reacs, on='id', how='right')

        # Get amount of missing reactions that have a MetaNetX ID
        statistics['mapped2MNX'] = len(biocyc_reacs.loc[biocyc_reacs['MetaNetX'].notna(), 'id'].unique().tolist())

        # @TODO
        # Split biocyc_reacs into MNX part & non-mapped part

        # @TODO
        # Step 2: Clean-up
        # ----------------
        # Move BioCyc ID to references & replace id with MetaNetX ID
    
        # Replace equation with mnx_equation if BioCyc could be mapped to MetaNetX

        return mnx_reacs, unmapped_reacs
    
    def _map_to_bigg(biocyc_reacs: pd.DataFrame) -> pd.DataFrame:
        """Maps Biocyc IDs in a table to the corresponding BiGG IDs & 
        Adds information obtained via the BiGG IDs to the resulting table

        Args:
            - biocyc_reacs (pd.DataFrame): 
                Table containing BioCyc IDs

        Returns:
            pd.DataFrame: 
                Table with BioCyc IDs mapped to BiGG IDs & more information 
                obtained from BiGG
        """

        # Step 1: Mapping & Obtain information
        # ------------------------------------
        # Mapping to BiGG + addition of more information on BiGG Reactions
        # Get BiGG BioCyc
        bigg2biocyc_reacs = load_a_table_from_database('SELECT id, BioCyc, name, reaction_string from bigg_reactions')

        # Change column names to fit input table
        bigg2biocyc_reacs.rename({'id': 'BiGG', 'BioCyc': 'id'}, axis=1, inplace=True)

        # Crop table to contain BiGG IDs set per BioCyc ID
        bigg_as_list = bigg2biocyc_reacs.groupby('id')['BiGG'].apply(set).reset_index(name='BiGG')
        bigg2biocyc_reacs.drop('BiGG', axis=1, inplace=True)
        bigg2biocyc_reacs = bigg_as_list.merge(bigg2biocyc_reacs, on='id')

        # Drop duplicates to get the unique BioCyc IDs
        bigg2biocyc_reacs.drop_duplicates(subset='id', inplace=True)

        # @TODO
        # Remove BiGG IDs for reactions in wrong compartments

        # Merge bigg2biocyc_reacs with biocyc_reacs to get new information
        biocyc_reacs = bigg2biocyc_reacs.merge(biocyc_reacs, on='id')

        # Get amount of missing reactions that have a BiGG ID
        statistics['mapped2BiGG'] = len(biocyc_reacs.loc[biocyc_reacs['BiGG'].notna(), 'id'].unique().tolist())

        # @TODO
        # Step 2: Clean-up
        # ----------------
        # Move BioCyc ID to references & replace id with BiGG ID
    
        # Replace equation with reaction_string if BioCyc could be mapped to BiGG

        # @DEPRECATE: Transfer this to model filling!
        # Add reactants & products dictionary with stoichiometric values to the reactions table
        # biocyc_reacs = add_stoichiometric_values_to_reacs_from_bigg(biocyc_reacs)

        return biocyc_reacs

    mapped2MNX = _map_to_mnx(biocyc_reacs)
    # mapped2BiGG = _map_to_bigg(mapped2MNX[mapped2MNX['via'] == 'BioCyc'])

    # @NOTE
    # @TODO
    missing_bc_reacs = mapped2MNX

    return missing_bc_reacs, statistics


# Mapping of EC numbers
# ---------------------

def map_ec_to_reac(table:pd.DataFrame, 
                   use_MNX:bool=True, use_BiGG:bool=True, 
                   use_KEGG:bool=True) -> pd.DataFrame:
    """Based on a table of NCBI protein IDs and EC numbers, 
    map them to reactions via different databases 
    (if a mapping is possible).
    
    input table should have format: 
        ec-code | ncbiprotein
        
    output table has the format:
        ec-code | ncbiprotein | id | equation | reference | is_transport | via

    Args:
        - table (pd.DataFrame): 
            The input table.
        - use_MNX (bool, optional): 
            Try mapping using the MetaNetX database. 
            Defaults to True.
        - use_BiGG (bool, optional): 
            Try mapping using the BiGG database. 
            Defaults to True.
        - use_KEGG (bool, optional): 
            Try mapping using the KEGG database. 
            Defaults to True.

    Returns:
        pd.DataFrame: 
            The extended table.
    """
    
    def _map_ec_to_reac_mnx(unmapped_reacs:pd.DataFrame) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the MetaNetX database.

        Args:
            - unmapped_reacs (pd.DataFrame): 
                The input table (ec-code and ncbiprotein columns)

        Returns:
            pd.DataFrame: 
                The extended table
        """
    
        # input: pd.DataFrame with at least the ec-code column
        # load MNX reac prop table
        mnx_reac_prop = load_a_table_from_database('mnx_reac_prop',False)
        # convert table into one EC-number per row
        mnx_reac_prop.drop('is_balanced', inplace=True, axis=1)
        mnx_reac_prop['ec-code'] = mnx_reac_prop['ec-code'].apply(lambda x: x.split(';') if isinstance(x,str) else None)
        # exclude entries without EC-number
        mnx_reac_prop = mnx_reac_prop.explode('ec-code').dropna(subset='ec-code')
        # merge with unmapped reactions
        reacs_mapped = unmapped_reacs.merge(mnx_reac_prop, on='ec-code', how='left')
        reacs_mapped.rename({'mnx_equation':'equation'}, inplace=True, axis=1)
        reacs_mapped['via'] = reacs_mapped['id'].apply(lambda x: 'MetaNetX' if x else None)
        
        return reacs_mapped      
    

    def _map_ec_to_reac_bigg(unmapped_reacs:pd.DataFrame) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the BiGG database.

        Args:
            - unmapped_reacs (pd.DataFrame): 
                The input table (ec-code and ncbiprotein columns)

        Returns:
            pd.DataFrame: 
                The extended table
        """
        
        # load BiGG reaction namespace
        bigg_reacs = load_a_table_from_database('bigg_reactions',False)
        bigg_reacs.dropna(subset='EC Number', inplace=True)
        bigg_reacs = bigg_reacs[['id','reaction_string','EC Number']].rename({'reaction_string':'equation','EC Number':'ec-code'}, inplace=False, axis=1)
        bigg_reacs['ec-code'] = bigg_reacs['ec-code'].apply(lambda x: x.split(',') if isinstance(x,str) else None)
        bigg_reacs = bigg_reacs.explode('ec-code')
        # merge with unmapped reactions
        bigg_mapping = unmapped_reacs.merge(bigg_reacs, on=['ec-code'], how='left')
        bigg_mapping.mask(bigg_mapping.isna(), other=None, inplace=True)
        # make conform to format
        bigg_mapping['reference'] = None
        bigg_mapping['is_transport'] = None
        bigg_mapping['via'] = bigg_mapping['id'].apply(lambda x: 'BiGG' if x else None)
        
        return bigg_mapping


    def _map_ec_to_reac_kegg(unmapped_reacs:pd.DataFrame) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the KEGG database.

        Args:
            - unmapped_reacs (pd.DataFrame): 
                The input table (ec-code and ncbiprotein columns)

        Returns:
            pd.DataFrame: 
                The extended table
        """
        
        # get KEGG EC number information
        eccodes = unmapped_reacs['ec-code']
        kegg_mapped = pd.DataFrame.from_dict(list(eccodes.progress_apply(parse_KEGG_ec)))
        kegg_mapped = unmapped_reacs.merge(kegg_mapped, on='ec-code')
        
        kegg_mapped['is_transport'] = None
        kegg_mapped['via'] = kegg_mapped['id'].apply(lambda x: 'KEGG' if x else None)
        kegg_mapped.explode(column='id')
        
        return kegg_mapped


    # input table should have format: 
    #   ec-code | ncbiprotein
    #   one EC number per row, list of ncbiprotein per row allowed
    if len(table.columns) != 2 or 'ec-code' not in table.columns or 'ncbiprotein' not in table.columns:
        raise ValueError('Wrong table format. Cannot map EC to reaction.')
    
    # map to MetaNetX
    if use_MNX:
        table = _map_ec_to_reac_mnx(table)
    
    # map to BiGG
    if use_BiGG:
        if 'id' in table.columns:
            to_map = table[table['id'].isna()][['ec-code','ncbiprotein']]    
            if len(to_map) > 0:        
                to_map = _map_ec_to_reac_bigg(to_map) 
                table = pd.concat([to_map,table[~table['id'].isna()]])
        else:
            table = _map_ec_to_reac_bigg(table)
            
    
    # map to KEGG
    if use_KEGG: 
        if 'id' in table.columns:
            to_map = table[table['id'].isna()][['ec-code','ncbiprotein']] 
            if len(to_map) > 0:
                to_map = _map_ec_to_reac_kegg(to_map) 
                table = pd.concat([to_map,table[~table['id'].isna()]])
        else:
            table = _map_ec_to_reac_kegg(table)
    
    # explode
    table = table.explode('id', ignore_index=True).explode('id', ignore_index=True)
        
    # output:   ec-code ncbiprotein	id	equation	reference	is_transport	via
    return table 


############################################################################
# classes
############################################################################

# -------------
# abtract class
# -------------

class GapFiller(ABC):
    """Abstract base class for the gap filling. 
    
    Already includes functions for the "filling" part of the gap-filling approach
    and some helper functions. Each subclass needs an implementation of `find_missing_genes` 
    and `find_missing_reactions` to determine the entities, that are missing in the model.

    Attributes:
        - full_gene_list (list): 
            List of all the genes. @TODO see code
        - geneid_type (str):
            What type of gene ID the model contains. 
            Defaults to 'ncbi'.
        - _statistics (dict):
            Dictionary of statistical information of the gap-filling run. Includes e.g.
            the number of added genes and reactions. 
    """

    def __init__(self) -> None:
        
        # data
        self.full_gene_list = None  # @TODO really a good idea? GeneGapFiller does not need it at all
        self.missing_genes = None       # missing genes, that have not yet been sorted into any category
        self.missing_reactions = None   # missing reacs, that have not yet been sorted into any category
        
        # general information
        self.geneid_type = 'ncbi' # @TODO more options?
        
        # collect stats & Co, can be extended by subclasses
        self._statistics = {  
                'genes':{
                    'missing (before)': 0,
                    'duplicates': 0,
                    'added': 0,
                    'missing (after)': 0
                },
                'reactions':{
                    'added (total)': 0,
                    'failed to build': 0
                },
                'metabolites':{
                    
                }
            }
        self.manual_curation = dict()
    
    # abstract methods
    # ----------------

    @abstractmethod
    # TODO write output format and what this functions needs to do
    def find_missing_genes(self, model):
        pass
    
    @abstractmethod
    # TODO write output format and what this functions needs to do
    def find_missing_reacs(self,model):
        pass
    
    # finding the gaps
    # ----------------

    def _find_reac_in_model(self, model: cobra.Model, eccode:str, id:str, 
               idtype:Literal['MetaNetX','KEGG','BiGG', 'BioCyc'], 
               include_ec_match:bool=False) -> Union[None, list]:
        """Helper function of :py:class:`~refinegems.classes.gapfill.GapFiller`.
        Search the model for an ID (and optionally EC number), to determine, if the 
        reaction is in the model.

        Args:
            - model (cobra.Model): 
                The model loaded with COBRApy.
            - eccode (str): 
                The EC number in the format: X.X.X.X
            - id (str): 
                The ID to search for.
            - idtype (Literal['MetaNetX','KEGG','BiGG', 'BioCyc']): 
                Name of the database the ID belongs to.
            - include_ec_match (bool, optional): 
                Option to include a match if only the EC number matches. 
                Defaults to False.

        Returns:
            Case one or more match found:

                list:

                    List of the ID of the reactions in the model, that match 
                    the query.
                    
            Case no match:
            
                None:
                
                    Nothing found.
        """
        
        # @TODO Ensure that user has requested BioCyc identifiers.org version? 
        # -> Could be done with polish_annotations
        MAPPING = {
            'MetaNetX':'metanetx.reaction', 
            'KEGG':'kegg.reaction',
            'BiGG':'bigg.reaction',
            'BioCyc': 'metacyc.reaction'
            }

        found = []
        for r in model.reactions:
            if MAPPING[idtype] in r.annotation.keys():
                if (isinstance(r.annotation[MAPPING[idtype]],list) and 
                    id in r.annotation[MAPPING[idtype]]):
                    found.append(r.id)
                elif (isinstance(r.annotation[MAPPING[idtype]],str) and 
                      id == r.annotation[MAPPING[idtype]]):
                    found.append(r.id)
            if include_ec_match and eccode and 'ec-code' in r.annotation.keys():
                if (isinstance(r.annotation['ec-code'],list) and 
                    eccode in r.annotation['ec-code']):
                    found.append(r.id)
                elif (isinstance(r.annotation['ec-code'],str) and 
                      eccode == r.annotation['ec-code']):
                    found.append(r.id)
            
        found = list(set(found))
        
        if len(found) > 0:
            return found
        
        return None
    
    
    # actual "Filling" part
    # ---------------------
    
    # @TODO logging? 
    def add_genes_from_table(self,model:libModel, gene_table:pd.DataFrame) -> None:
        """Create new GeneProduct for a table of genes in the format:
        
        | ncbiprotein | locus_tag | UniProt | ... |
        
        The dots symbolise additional columns, that can be passed to the function, 
        but will not be used by it. The other columns, except UniProt, are required.

        Args:
            - model (libModel): 
                _description_
            - gene_table (pd.DataFrame): 
                The table with the genes to add. At least needs the columns
                '','' and 'ec-code'. 
                @TODO: Hier fehlt was, oder? Passt auch 
                nicht ganz zu den Spaltennamen oben...
        """
    
        # ncbiprotein | locus_tag | ...
        # work on a copy to ensure input stays the same
        gene_table = gene_table.copy()
        # gene_table.drop(columns=['ec-code'],inplace=True)
        
        # create gps from the table and add them to the model
        if 'UniProt' in gene_table.columns:
            for idx,x in tqdm(gene_table.iterrows(), 
                              total=len(gene_table),
                              desc='Adding genes to model'):
                create_gp(model, x['ncbiprotein'], 
                        locus_tag=x['locus_tag'],
                        uniprot=(x['UniProt'],True))
        else:
            for idx,x in tqdm(gene_table.iterrows(), 
                              total=len(gene_table),
                              desc='Adding genes to model'):
                create_gp(model, x['ncbiprotein'], 
                        locus_tag=x['locus_tag'])
                

    # @TODO seems very rigid, better ways to find the ids?
    def add_gene_reac_associations_from_table(self,model:libModel,
                                              reac_table:pd.DataFrame) -> None:
        """Using a table with at least the columns 'ncbiprotein' 
        (containing e.g. NCBI protein identifier (lists), should be gene IDs in the model)
        and 'add_to_GPR' (containing reactions identifier (lists)), add the gene IDs to the 
        GPRs of the corresponding reactions.

        Args:
            - model (libModel): 
                The model loaded with libSBML.
            - reac_table (pd.DataFrame): 
                The table containing at least the columns 'ncbiprotein' (gene IDs) and
                'add_to_GPR' (reaction IDs)
        """
        
        model_gene_ids = [_.getId() for _ in model.getPlugin(0).getListOfGeneProducts()]
        
        # get each unique ncbiprotein vs reaction mapping
        reac_table = reac_table[['ncbiprotein','add_to_GPR']]
        reac_table = reac_table.explode('ncbiprotein').explode('add_to_GPR')
        reac_table.drop_duplicates(inplace=True)
        
        # add the genes to the corresponding GPRs
        for idx,row in reac_table.iterrows():
            # check, if G_+ncbiprotein in model
            # if yes, add gpr
            geneid = 'G_'+row['ncbiprotein'].replace('.','_')
            reacid = 'R_'+row['add_to_GPR']
            if geneid in model_gene_ids:
                create_gpr(model.getReaction(reacid),geneid)
            # else, print warning
            else:
                mes = f'Cannot find {geneid} in model. Should be added to {reacid}'
                warnings.warn(mes,UserWarning)
    

    # @TEST BioCyc
    # @TODO logging, save stuff for manual curation etc. -> started a bit
    # @TEST - somewhat seems to work - for now
    def add_reactions_from_table(self, model:cobra.Model,
                                 missing_reac_table:pd.DataFrame,
                                 formula_check:Literal['none','existence','wildcard','strict']='existence',
                                 exclude_dna:bool=True,
                                 exclude_rna:bool=True,
                                 idprefix:str='refineGEMs',
                                 namespace:Literal['BiGG']='BiGG') -> pd.DataFrame:
        """Helper function to add reactions to a model from the missing_reactions table 
        (output of the chosen implementation of :py:meth:`~refinegems.classes.gapfill.GapFiller.find_missing_reacs`)

        Args:
            - model (cobra.Model): 
                The model, loaded with COBRpy.
            - missing_reac_table (pd.DataFrame): 
                The missing reactions table.  
            - formula_check (Literal['none','existence','wildcard','strict'], optional):
                Param for checking metabolite formula before adding them to the model.
                For more information, refer to :py:func:`~refinegems.utility.entities.isreaction_complete`
                Defaults to 'existence'.
            - exclude_dna (bool, optional):
                Option to exclude reaction containing 'DNA' from being added to the model.
                Defaults to True.
            - exclude_rna (bool, optional):
                Option to exclude reaction containing 'RNA' from being added to the model.
                Defaults to True.
            - idprefix (str, optional): 
                A prefix to use, if pseudo-IDs need to be created. 
                Defaults to 'refineGEMs'.
            - namespace (Literal['BiGG'], optional): 
                Namespace to use for the reactions and metabolites 
                (and the model). Defaults to 'BiGG'.
                
        Raises:
            - TypeError: Unknown return type for reac param. Please contact the developers.

        Returns:
            pd.DataFrame: 
                Table containing the information about which genes can now be added 
                to reactions (use for GPR curation).
        """
        
        # reconstruct reactions
        # ---------------------
        for idx,row in tqdm(missing_reac_table.iterrows(), 
                            desc='Trying to add missing reacs',
                            total=missing_reac_table.shape[0]):
            # build reaction
            reac = None
            match row['via']:
                # MetaNetX
                case 'MetaNetX':
                    reac = build_reaction_mnx(model,row['id'],
                                            reac_str=row['equation'],
                                            references={'ec-code':[row['ec-code']]},
                                            idprefix=idprefix,
                                            namespace=namespace)      
                # KEGG
                case 'KEGG':
                    refs = row['references']
                    refs['ec-code'] = row['ec-code']
                    reac = build_reaction_kegg(model,row['id'],
                                            reac_str=row['equation'],
                                            references=refs,
                                            idprefix=idprefix,
                                            namespace=namespace)
                # BiGG
                case 'BiGG':
                    reac = build_reaction_bigg(model,row['id'],
                                            references={'ec-code':[row['ec-code']]},
                                            idprefix=idprefix,
                                            namespace=namespace)
                    
                # Unknown database
                case _:
                    mes = f'''Unknown database name for reaction reconstruction: {row["via"]}\n
                    Reaction will not be reconstructed.'''
                    warnings.warn(mes,UserWarning)
            
            # check output of reconstruction
            # ------------------------------
            # case 1: reconstruction was not possible
            if not reac:
                pass # nothing to do here
            # case 2: reaction(s) found in model
            elif isinstance(reac,list):
                # add found names to the add_to_GPR column of the table
                current_gpr = missing_reac_table.loc[idx,'add_to_GPR']
                if not current_gpr:
                    missing_reac_table.at[idx,'add_to_GPR'] = reac
                else:
                    missing_reac_table.at[idx,'add_to_GPR'] = list(set(reac + current_gpr))
            # case 3: new reaction was generated
            elif isinstance(reac,cobra.Reaction):
                # validate reaction
                if isreaction_complete(reac, formula_check=formula_check,
                                       exclude_dna=exclude_dna,
                                       exclude_rna=exclude_rna):
                    # add reaction to model (if validation succesful)
                    model.add_reactions([reac])
                    self._statistics['reactions']['added (total)'] = self._statistics['reactions']['added (total)'] + 1
                    # add reaction ID to table under add_to_GPR
                    current_gpr = missing_reac_table.loc[idx,'add_to_GPR']
                    if not current_gpr:
                        missing_reac_table.at[idx,'add_to_GPR'] = [reac.id]
                    else:
                        current_gpr.append(reac.id)
                        missing_reac_table.at[idx,'add_to_GPR'] = list(set(current_gpr))
            # case 4: should never occur
            else:
                mes = f'Unknown return type for reac param. Please contact the developers.'
                raise TypeError(mes)
            
        # save reactions, that could not be recontructed, for manual curation
        manual_curation_reacs = missing_reac_table[missing_reac_table['add_to_GPR'].isnull()]
        self.manual_curation['reaction, building failed'] = manual_curation_reacs
        self._statistics['reactions']['failed to build'] = len(manual_curation_reacs)
        
        # return the updated table with successfully reconstructed reaction ids 
        # to enable adding the genes
        missing_gprs = missing_reac_table[~missing_reac_table['add_to_GPR'].isnull()]
        return missing_gprs


    # @TODO : logging
    # @TODO : save stuff for report / manual curation
    # @DISCUSSION: Check for futile cycles after each addition of a new reaction?
    # @TEST - only tested in debugging mode with GeneGapFiller
    def fill_model(self, model:Union[cobra.Model,libModel], 
                   **kwargs) -> libModel:
        """Based on a table of missing genes and missing reactions, 
        fill the gaps in a model as good as possible automatically.
        
        .. note::
        
            This model rewrites and reloads the input model. Only the returned model
            has all the edits.

        Args:
            - model (Union[cobra.Model,libModel]): 
                The model, either a libSBML or COBRApy model entity.
            - kwargs:
                Additional parameters to be passed to 
                :py:meth:`~refinegems.classes.gapfill.GapFiller.add_reactions_from_table`.

        Raises:
            - TypeError: Unknown type of model.

        Returns:
            libModel: 
                The gap-filled model.
        """
        # filter out duplicates genes to avoid duplicates IDs in the model
        # @TODO: Better idea to add duplicate genes
        if len(self.missing_genes) != len(self.missing_genes['ncbiprotein'].unique()):
            self.manual_curation['duplicate genes (not added)'] = self.missing_genes[self.missing_genes.duplicated(subset=['ncbiprotein'])]
            self._statistics['genes']['duplicates'] = self._statistics['genes']['duplicates'] + len(self.manual_curation['duplicate genes (not added)'])
            self.missing_genes = self.missing_genes[~self.missing_genes.duplicated(subset=['ncbiprotein'])]
    
        
        # load the correct type of model for the first step
        # @TODO probably not working under Windows
        match model:
            case cobra.Model():
                with NamedTemporaryFile(suffix='.xml') as tmp:
                    write_model_to_file(model,tmp.name)
                    model = load_model(tmp.name,'libsbml')
            case libModel():
                pass
            case _:
                mes = f'Unknown type of model: {type(model)}'
                raise TypeError(mes)
        
        # Step 1: Add genes to model whose reactions are already in it
        # -------------------------------------------------------------
        # filter the respective genes and reactions
        reacs_in_model = self.missing_reacs[~(self.missing_reacs['add_to_GPR'].isnull())]
        ncbiprot_with_reacs_in_model = [*chain(*list(reacs_in_model['ncbiprotein']))]
        genes_with_reacs_in_model = self.missing_genes[self.missing_genes['ncbiprotein'].isin(ncbiprot_with_reacs_in_model)]
        self._statistics['genes']['added'] = self._statistics['genes']['added'] + len(genes_with_reacs_in_model)
        if len(genes_with_reacs_in_model) > 0:
            # add genes as gene products to model
            self.add_genes_from_table(model, genes_with_reacs_in_model)
        
            # extend gene production rules 
            self.add_gene_reac_associations_from_table(model,reacs_in_model)
            
            # what remains:
            self.missing_reacs = self.missing_reacs[self.missing_reacs['add_to_GPR'].isnull()]
            self.missing_genes = self.missing_genes[~(self.missing_genes['ncbiprotein'].isin(ncbiprot_with_reacs_in_model))]
        
        
        # Step 2: Add reactions to model, if reconstruction successful
        # ------------------------------------------------------------

        if len(self.missing_reacs) > 0:
            
            # re-load model with cobrapy
            with NamedTemporaryFile(suffix='.xml') as tmp:
                write_model_to_file(model,tmp.name)
                cobramodel = load_model(tmp.name,'cobra')
            
            # .......................
            # @DEBUG
            # if len(self.missing_reacs) > 10:
            #     self.missing_reacs = self.missing_reacs.sample(10)
            #     print('fill_model: Running in debugging mode')
            # .......................
                
            # add reactions to model  
            missing_gprs = self.add_reactions_from_table(cobramodel,self.missing_reacs,**kwargs)

        # Step 3: Add GPRs + genes for the newly curated reactions 
        # --------------------------------------------------------
        
        # re-load model with libsbml
        # @TODO does not seem to work as expected under Windows
        with NamedTemporaryFile(suffix='.xml') as tmp:
            write_model_to_file(cobramodel,tmp.name)
            model = load_model(tmp.name,'libsbml')
            
        if len(missing_gprs) > 0:
            # filter for genes for GPRs but not yet in model
            ncbiprot_with_reacs_in_model = [*chain(*list(missing_gprs['ncbiprotein']))]
            genes_with_reacs_in_model = self.missing_genes[self.missing_genes['ncbiprotein'].isin(ncbiprot_with_reacs_in_model)]
            self._statistics['genes']['added'] = self._statistics['genes']['added'] + len(genes_with_reacs_in_model)
            if len(genes_with_reacs_in_model) > 0:
                # add genes as gene products to model
                self.add_genes_from_table(model, genes_with_reacs_in_model)
                # extend gene production rules 
                self.add_gene_reac_associations_from_table(model,reacs_in_model)
        
        # collect stats and stuff for manual curation
        self.missing_genes = self.missing_genes[~(self.missing_genes['ncbiprotein'].isin(ncbiprot_with_reacs_in_model))]
        self.manual_curation['missing genes (after gap filling)'] = self.missing_genes
        self._statistics['genes']['missing (after)'] = len(self.missing_genes)
        
        return model
        

    
    # reporting -> new class in reports?
    # ---------

    def calculate_stats(self):
        pass

    def report(self):
        pass

# --------------------
# Gapfilling with KEGG
# --------------------

# @TEST
class KEGGapFiller(GapFiller):
    """Based on a KEGG organism ID (corresponding to the organism of the model),
    find missing genes in the model and map them to reactions to try and fill the gaps
    found with the KEGG database. 
    
    .. note::
    
        Due to the KEGG REST API this is relatively slow.

    Attributes:
    
    - GapFiller Attributes:
        all attributes of the parent class, :py:class:`~refinegems.classes.gapfill.GapFiller`
    - organismid (str, required):
        Abbreviation of the organism in the KEGG database.
    
    """

    def __init__(self, organismid) -> None:
        super().__init__()
        self.organismid = organismid
        
        
    # @TODO: parallelising
    # @TODO: logging
    def find_missing_genes(self, model:libModel):
        """Get the missing genes in model in comparison to the KEGG entry of the 
        organism. Saves a table containing the missing genes (format: <organism id>:<locus>)
        to the attribute missing_genes.

        Args:
            - model (libModel): 
                The model loaded with libSBML.
                
        """
    
        # Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis/entities --- Modified
        def get_model_genes(model: libModel) -> pd.DataFrame:
            """Extracts KEGG Genes from given model

            Args:
                - model (model-libsbml): 
                    Model loaded with libSBML

            Returns:
                pd.DataFrame: 
                    Table with all KEGG Genes in the model
            """
            genes_in_model = []
            for gene in model.getPlugin(0).getListOfGeneProducts():
                cv_terms = gene.getCVTerms()
                if cv_terms:
                    for cv_term in cv_terms:
                        for idx in range(cv_term.getNumResources()):
                            uri = cv_term.getResourceURI(idx)
                            if 'kegg.genes' in uri: 
                                genes_in_model.append(re.split('kegg.genes:|kegg.genes/',uri)[1]) # work with old/new pattern

            return pd.DataFrame(genes_in_model, columns=['orgid:locus'])
        
        
        # Step 1: get genes from model
        # ----------------------------
        genes_in_model = get_model_genes(model)
        
        # Step 2: get genes of organism from KEGG
        # ---------------------------------------
        gene_KEGG_list = KEGG().list(self.organismid)
        gene_KEGG_table = pd.read_table(io.StringIO(gene_KEGG_list), header=None)
        gene_KEGG_table.columns = ['orgid:locus','CDS','position','protein']
        gene_KEGG_table = gene_KEGG_table[['orgid:locus']]
        
        # Step 3: KEGG vs. model genes -> get missing genes for model
        # ----------------------------
        genes_not_in_model = gene_KEGG_table[~gene_KEGG_table['orgid:locus'].isin(genes_in_model['orgid:locus'])]
        
        # Step 4: extract locus tag
        # -------------------------
        genes_not_in_model['locus_tag'] = genes_not_in_model['orgid:locus'].str.split(':').str[1]
        
        # Step 5: map to EC via KEGG
        # --------------------------
        # @DEBUG .......................
        # genes_not_in_model = genes_not_in_model.iloc[330:350,:]
        # print(UserWarning('Running in debugging mode.'))
        # ..............................
        geneKEGG_mapping = pd.DataFrame.from_dict(list(genes_not_in_model['orgid:locus'].progress_apply(parse_KEGG_gene)))
        genes_not_in_model = genes_not_in_model.merge(geneKEGG_mapping, how='left', on='orgid:locus')
        genes_not_in_model = genes_not_in_model.explode('ncbiprotein')
        
        # collect stats
        self._statistics['genes']['missing (before)'] = len(genes_not_in_model)
        
        self.missing_genes = genes_not_in_model 
    
    
    # @TODO : logging
    # @TODO : paralellising possibilities?
    # @TODO : progress bar
    # @TODO : self._statistics for reactions
    def find_missing_reacs(self,model:cobra.Model):
 
        # Step 1: filter missing gene list + extract ECs
        # ----------------------------------------------
        reac_options = self.missing_genes[['ec-code','ncbiprotein']]        # get relevant infos for reacs
        self.missing_reacs = reac_options[['ec-code','ncbiprotein']].dropna()    # drop nas
        self.manual_curation['no EC/ncbiprotein'] = reac_options.loc[~reac_options.index.isin(self.missing_reacs.index)]
        # check, if any automatic gapfilling is possible
        if len(self.missing_reacs) == 0:
            return None
        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = self.missing_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = self.missing_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        self.missing_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        self.missing_reacs = self.missing_reacs.groupby(self.missing_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
        
        # Step 2: map EC to reaction(s) if possible
        # -----------------------------------------
        # via MNX, BiGG, KEGG
        reacs_mapped = map_ec_to_reac(self.missing_reacs)
        
        # Step 3: clean and map to model reactions
        # ----------------------------------------
        # need manual curation
        self.manual_curation['reacs'] = reacs_mapped[reacs_mapped['id'].isnull()]
        # map to model reactions
        reacs_mapped = reacs_mapped[~reacs_mapped['id'].isnull()] 
        reacs_mapped['add_to_GPR'] = reacs_mapped.apply(lambda x: self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1)
                
        self.missing_reacs = reacs_mapped
    
    
# ----------------------
# Gapfilling with BioCyc
# ----------------------
# @TODO: Possibility to add KEGG Reaction & MetaNetX IDs to SmartTable 
#       -> So far these are empty for my test cases
#       -> Could be added in a future update
class BioCycGapFiller(GapFiller):
    """
    | @TODO: Write correct doc string
    | These three TXT files can be obtained through creating SmartTables in BioCyc 
    | and exporting these as 
    | Spreadsheets with the parameter FrameID (one & two) or Common name (three).
    | SmartTable three: Should contain 'Compound', 'Chemical Formula' and 'InChI-Key'.

    Attributes:
        - biocyc_gene_tbl_path (str, required): 
            Path to organism-specific SmartTable for genes from BioCyc;
            Should contain the columns: 'Accession-2' | 'Reactions of gene'
        - biocyc_reacs_tbl_path (str, required): 
            Path to organism-specific SmartTable for reactions from BioCyc;
            Should contain the columns: 
            'Reaction' | 'Object ID' | 'EC-Number' | 'Spontaneous?'
        - gff (str, required): 
            Path to organism-specific GFF file
        - missing_genes (pd.DataFrame):
            DataFrame containing the missing genes with the columns 
            'locus_tag' | 'id' | 'ncbiprotein' | 'name'
        - missing_reacs (pd.DataFrame):
            DataFrame containing the missing reactions with the columns 
            'id' | 'ncbiprotein' | 'equation' | 'ec-code' | 'via' | 'add_to_GPR'
    """
    
    def __init__(self, biocyc_gene_tbl_path: str, 
                 biocyc_reacs_tbl_path: str, gff:str) -> None:
        super().__init__()
        self.biocyc_gene_tbl = biocyc_gene_tbl_path
        # @TODO: Das ist eigtl self.full_gene_list von GapFiller!
        self.biocyc_rxn_tbl = biocyc_reacs_tbl_path
        self._gff = gff
        self.missing_genes = None
        self.missing_reacs = None

    @property
    def biocyc_gene_tbl(self):
        return self._biocyc_gene_tbl
    
    @biocyc_gene_tbl.setter
    # @DISCUSSION: Hier sollten wir noch diskutieren, ob Accession-2 oder Accession-1 
    # hard coded sein sollten oder nicht
    #        Dafür müssten wir nochmal ein paar entries in BioCyc dazu anschauen.
    # Locus tags in GenBank GFF file == BioCyc Accession-2 == Old locus tags in 
    # RefSeq GFF file
    # Locus tags in RefSeq GFF file == BioCyc Accession-1
    # Label in model == Locus tag from GenBank GFF file == BioCyc Accession-2
    def biocyc_gene_tbl(self, biocyc_gene_tbl_path: str):
        """Parses TSV file from BioCyc to retrieve 'Accession-2' & the 
        corresponding 'Reactions of gene'
    
        Args:
            - inpath (str): 
                Path to organism-specific SmartTable for genes from BioCyc;
                Should contain the columns: 'Accession-2' | 'Reactions of gene'
           
        Returns:
            pd.DataFrame: 
                Table containing only rows where a 'Reaction of gene' exists
        """
        # Read table
        self._biocyc_gene_tbl = pd.read_table(
            biocyc_gene_tbl_path, 
            usecols=['Accession-2', 'Reactions of gene'], 
            dtype=str
            )

        # Rename columns for further use
        self._biocyc_gene_tbl.rename(
            columns={
                'Accession-2': 'locus_tag', 'Reactions of gene': 'id'
                }, 
            inplace=True)

        # Turn empty strings into NaNs
        self._biocyc_gene_tbl.replace('', np.nan, inplace=True)

        # Drop only complete empty rows
        self._biocyc_gene_tbl.dropna(how='all', inplace=True)

        # Save not mappable genes
        self.manual_curation['BioCyc genes unmappable'] = self._biocyc_gene_tbl[self._biocyc_gene_tbl['id'].isna()]

        # Add amount of unmappable genes to statistics
        self._statistics['genes']['missing (unmappable)'] = len(
            self.manual_curation['BioCyc genes unmappable']['locus_tag'].unique().tolist()
            )

        # Remove all rows where 'id' NaNs
        self._biocyc_gene_tbl.dropna(subset='id', inplace=True)

    @property
    def biocyc_rxn_tbl(self):
        return self._biocyc_rxn_tbl

    @biocyc_rxn_tbl.setter
    def biocyc_rxn_tbl(self, biocyc_reacs_tbl_path: str) -> pd.DataFrame:
        """Parses TSV file from BioCyc to retrieve 'Reaction', 'Object ID', 
        'EC-Number' & 'Spontaneous?'

        Args:
            - biocyc_reacs_tbl_path (str):   
                Path to organism-specific SmartTable for reactions from BioCyc;
                Should contain the columns: 
                'Reaction' | 'Object ID' | 'EC-Number' | 'Spontaneous?'

        Returns:
            pd.DataFrame: 
                Table containing all BioCyc reactions from provided file
        """
        # Read table
        self._biocyc_rxn_tbl = pd.read_table(
            biocyc_reacs_tbl_path, 
            usecols=[
                'Reaction', 'Object ID', 'EC-Number', 'Spontaneous?'
                ],
            dtype=str
            )

        # Rename columns for further use
        self._biocyc_rxn_tbl.rename(
            columns={
                'Reaction': 'equation', 'Object ID': 'id',
                'EC-Number': 'ec-code', 'Spontaneous?': 'is_spontaneous'
                },
            inplace=True
            )

        # Turn empty strings into NaNs
        self._biocyc_rxn_tbl.replace('', np.nan, inplace=True)

        # Set entries in is_spontaneous to booleans &
        # specify empty entries in 'is_spontaneous' as False
        self._biocyc_rxn_tbl['is_spontaneous'].replace({'T': True, 'F': False}, inplace=True)
        self._biocyc_rxn_tbl['is_spontaneous'] = self._biocyc_rxn_tbl['is_spontaneous'].fillna(False)

    def find_missing_genes(self, model: libModel):
        """Retrieves the missing genes and reactions from the BioCyc table 
        according to the 'Accession-2' identifiers

        Args:
            - model (libModel): 
                Model loaded with libSBML
        """
        
        # Step 1: get genes from model
        # ----------------------------
        geneps_in_model = [
            _.getLabel() 
            for _ in model.getPlugin(0).getListOfGeneProducts()
            ]

        # Step 2: Get genes of organism from BioCyc
        # -----------------------------------------
        # See setter: BioCyc_gene_tbl
        
        # Step 3: BioCyc vs. model genes -> get missing genes for model
        # -------------------------------------------------------------
        self.missing_genes = self.biocyc_gene_tbl[
            ~self.biocyc_gene_tbl['locus_tag'].isin(geneps_in_model)
            ]

        # Step 4: Get ncbiprotein IDs
        # ---------------------------
        # Parse GFF file to obtain locus_tag2ncbiportein mapping for all CDS
        locus_tag2ncbiprotein_df = parse_gff_for_cds(
            self._gff,
            {
                'locus_tag': 'locus_tag', 
                'protein_id': 'ncbiprotein',
                'product': 'name'}
            )
        locus_tag2ncbiprotein_df = locus_tag2ncbiprotein_df.explode('ncbiprotein')
        locus_tag2ncbiprotein_df = locus_tag2ncbiprotein_df.explode('name')

        # Get the complete missing genes dataframe with the ncbiprotein IDs
        self.missing_genes = self.missing_genes.merge(
            locus_tag2ncbiprotein_df, on='locus_tag'
            )

        # Step 5: Get amount of missing genes from BioCyc for statistics
        # --------------------------------------------------------------
        self._statistics['genes']['missing (before)'] = len(
            self.missing_genes['locus_tag'].unique().tolist()
            )

    def find_missing_reacs(self, model: cobra.Model):
        """Retrieves the missing reactions with more information like the 
        equation, EC code, etc. according to the missing genes

        Args:
            - model (cobra.Model): 
                Model loaded with COBRApy
        """

        # Step 1: filter missing gene list + extract ECs
        # ----------------------------------------------
        # Drop locus tag column as not needed for here
        missing_genes = self.missing_genes.drop(['locus_tag', 'name'], axis=1)
        
        # Expand missing genes result table to merge with Biocyc reactions table
        missing_genes = pd.DataFrame(
           missing_genes['id'].str.split('//').tolist(), 
           index=missing_genes['ncbiprotein']
           ).stack()
        missing_genes = missing_genes.reset_index(
            [0, 'ncbiprotein']
            )
        missing_genes.columns = ['ncbiprotein', 'id']
        missing_genes['id'] = missing_genes['id'].str.strip()

        # Turn ncbiprotein column into lists of ncbiprotein IDs per reaction
        ncbiprotein_as_list = missing_genes.groupby('id')['ncbiprotein'].apply(list).reset_index(name='ncbiprotein')
        missing_genes.drop('ncbiprotein', axis=1, inplace=True)
        missing_genes = ncbiprotein_as_list.merge(missing_genes, on='id')
        missing_genes['ncbiprotein'] = missing_genes['ncbiprotein'].apply(
            lambda x: 
                x if not None in x else list(filter(None, x))
        )

        # Drop duplicates to get the unique BioCyc IDs
        missing_genes.drop_duplicates(subset='id', inplace=True)

        # Get missing reactions from missing genes
        self.missing_reacs = missing_genes.merge(
            self.biocyc_rxn_tbl, on='id'
            )

        # Turn ec-code entries with '//' into lists
        self.missing_reacs['ec-code'] = self.missing_reacs['ec-code'].str.split('\s*//\s*')
        
        # Add 'G_spontaneous' as gene product if marked as spontaneous &
        # drop is_spontaneous column
        self.missing_reacs['ncbiprotein'] = self.missing_reacs.apply(
            lambda x: 
                x['ncbiprotein'] if not x['is_spontaneous'] else x['ncbiprotein'].append('spontaneous'),
            axis=1
        )
        self.missing_reacs.drop('is_spontaneous', axis=1, inplace=True)

        # Step 2: Get amount of missing reactions from BioCyc for statistics
        # ------------------------------------------------------------------
        self._statistics['reactions']['missing (before)'] = len(
            self.missing_reacs['id'].unique().tolist()
            )

        # Step 3: Map BioCyc to model reactions & cleanup
        # -----------------------------------------------
        # Add column 'via'
        self.missing_reacs['via'] = 'BioCyc'

        # Filter reacs for already in model
        self.missing_reacs['add_to_GPR'] = self.missing_reacs.apply(
            lambda x: 
                self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1
            )

        # Add column 'references'
        self.missing_reacs['references'] = None

        # @TODO
        # Step 4: Map missing reactions without entries in column 'add_to_GPR' 
        #         to other databases to get a parsable reaction equation
        # --------------------------------------------------------------------
        # Map to MetaNetX, then to BiGG
        mapped_reacs = None # map_biocyc_to_reac(self.missing_reacs[self.missing_reacs['add_to_GPR'].isna()])

        
        # Filter reacs for already in model
        '''
        mapped_reacs['add_to_GPR'] = mapped_reacs.apply(
            lambda x: 
                self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1
            )
        '''

        # Merge self.missing_reacs with the mapped_reacs
        # pd.concat(self.missing_reacs[self.missing_reacs['add_to_GPR'].isna()], mapped_reacs)

        # Step 5: Get results
        # -------------------
        # Split missing reactios based on entries in 'via' & 'add_to_GPR'
        mask = (self.missing_reacs['via'] == 'BioCyc') | (self.missing_reacs['add_to_GPR'].isna())
        
        # DataFrame with unmappable BioCyc IDs & No entries in 'add_to_GPR'
        self.manual_curation['BioCyc reactions unmappable'] = self.missing_reacs[mask]

        # DataFrame with either mapped BioCyc IDs or Entries in 'add_to_GPR'
        self.missing_reacs = self.missing_reacs[~mask]
    
# ----------------
# Gapfilling no DB
# ----------------
# @NOTE: Ideas from Gwendolyn O. Döbel for lab strains
# Get all  possible genes by filtering .gff according to 'bio_type=protein_coding' & 'product=hypothetical protein'
# Compare the list of genes with the ones already in the model & add all missing genes
# Before adding to model check if for all genes that are missing for IMITSC147 identifiers exist
# -> Create tables mapping locus tag to old ID, locus tag to new ID & merge 
# -> Specify user input locus_tag start from NCBI PGAP
# # Skeleton for functions that could be used for a lab strain/organism which is in no database contained
# def get_genes_from_gff():
#     pass
# 
# 
# def get_related_metabs_reactions_blast():
#     pass
# 
# 
# def gff_gene_comp():
#     pass
# 
#

# ---------------------------------
# GapFilling with GFF and Swissprot
# ---------------------------------

# @TEST
class GeneGapFiller(GapFiller):
    """Find gaps in the model using the GFF file of the underlying genome 
    and the Swissprot database and optionally NCBI.
    
    This gap filling approach tries to identify missing genes from the GFF file
    and uses DIAMOND to run a blastp search for homologs against the Swissprot database.
    
    .. hint::
        Files required for this approach can be downloaded with :py:func:`~refinegems.utility.set_up.download_url`
    
    """
    
    GFF_COLS = {'locus_tag':'locus_tag', 
                    'eC_number':'ec-code', 
                    'protein_id':'ncbiprotein'} # :meta: 
    
    def __init__(self) -> None:
        super().__init__()
        
    def find_missing_genes(self,gffpath:Union[str|Path],model:libModel):
    
        # get all CDS from gff
        all_genes = parse_gff_for_cds(gffpath,self.GFF_COLS)
        # get all genes from model by locus tag
        model_locustags = [g.getLabel() for g in model.getPlugin(0).getListOfGeneProducts()]
        # filter
        self.missing_genes = all_genes.loc[~all_genes['locus_tag'].isin(model_locustags)]
        # formatting
        for col in self.GFF_COLS.values():
            if col not in self.missing_genes.columns:
                self.missing_genes[col] = None

        # collect stats
        self._statistics['genes']['missing (before)'] = len(self.missing_genes)
                
        # save genes with no locus tag for manual curation
        self.manual_curation['gff no locus tag'] = self.missing_genes[self.missing_genes['locus_tag'].isna()]['ncbiprotein']
        self._statistics['genes']['no locus tag'] = len(self.manual_curation['gff no locus tag'])
        
        # formatting
        # ncbiprotein | locus_tag | ec-code
        self.missing_genes =  self.missing_genes[~self.missing_genes['locus_tag'].isna()]
        self.missing_genes = self.missing_genes.explode('ncbiprotein')
        
    
    def find_missing_reacs(self, model:cobra.Model,  
                          # prefix for pseudo ncbiprotein ids
                          prefix:str='refinegems',
                          # NCBI params
                          mail:str=None, 
                          check_NCBI:bool=False,
                          # SwissProt
                          fasta:str=None, 
                          dmnd_db:str=None, 
                          swissprot_map:str=None,
                          **kwargs) -> tuple:
        
        # Case 1:  no EC
        # --------------
        case_1 = self.missing_genes[self.missing_genes['ec-code'].isna()]
        not_case_1 = self.missing_genes[~self.missing_genes['ec-code'].isna()]
        if len(case_1) > 0:
            
            # Option 1: BLAST against SwissProt
            # +++++++++++++++++++++++++++++++++    
            # -> BLAST (DIAMOND) against SwissProt to get EC/BRENDA 
            if fasta and dmnd_db and swissprot_map:
                case_1_mapped = get_ec_via_swissprot(fasta,dmnd_db,
                                            case_1,
                                            swissprot_map,
                                            **kwargs) # further optional params for the mapping
                case_1.drop('ec-code', inplace=True, axis=1)
                case_1 = case_1.merge(case_1_mapped, on='locus_tag', how='left')
                not_case_1['UniProt'] = None
                
            # Option 2: Use ML to predict EC
            # ++++++++++++++++++++++++++++++
            # @TODO
            # -> sth like DeepECTransformer (tools either not good, 
            #    not installable or no license)
            # -> use ECRECer Web service output as input 
            #    (whole protein fasta -> wait for Gwendolyn's results) -> does not work / not return
            # -> use CLEAN webservice
            #    same problem as above with the web tool

        self.manual_curation['no ncbiprotein, no EC'] = case_1[case_1['ncbiprotein'].isna() & case_1['ec-code'].isna()] 
        
        mapped_reacs = pd.concat([case_1[~(case_1['ncbiprotein'].isna() & case_1['ec-code'].isna())],not_case_1])
        self._statistics['reactions']['no NCBI, no EC'] = len(self.manual_curation['no ncbiprotein, no EC'])

        # convert NaNs to None
        mapped_reacs.mask(mapped_reacs.isna(), other=None, inplace=True)

        # Case 2: still no EC but ncbiprotein
        # -----------------------------------
        #       -> access ncbi for ec (optional) 
        # @DEBUG .......................
        # mapped_reacs = mapped_reacs.iloc[300:350,:]
        # print(UserWarning('Running in debugging mode.'))
        # ..............................
        if check_NCBI and mail:
            mapped_reacs['ec-code'] = mapped_reacs.progress_apply(lambda x: get_ec_from_ncbi(mail,x['ncbiprotein']) if not x['ec-code'] and not x['ncbiprotein'].isna() else x['ec-code'], axis=1)
        
        # save entries with no EC for manual curation
        self.manual_curation['no EC'] = mapped_reacs[mapped_reacs['ec-code'].isna()]
        self._statistics['reactions']['NCBI, no EC'] = len(self.manual_curation['no EC'])
        mapped_reacs = mapped_reacs[~mapped_reacs['ec-code'].isna()]
        
        # check, if any automatic gapfilling is still possible
        if len(mapped_reacs) == 0:
            return None
        
        # create pseudoids for entries with no ncbiprotein id
        mapped_reacs['ncbiprotein'] = mapped_reacs.apply(lambda x: f'{prefix}_{x["locus_tag"]}' if not x['ncbiprotein'] else x['ncbiprotein'], axis=1)

        # Case 3: EC found
        # ----------------
        
        # update the gene information
        updated_missing_genes = mapped_reacs.copy()
        
        # reformat missing reacs 
        mapped_reacs.drop(['UniProt','locus_tag'], inplace=True, axis=1)
        
        # transform table into EC-number vs. list of NCBI protein IDs
        # @TODO make a func out of this - occurs on multiple occasions
        eccode = mapped_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = mapped_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        mapped_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        mapped_reacs = mapped_reacs.groupby(mapped_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
        
        # map EC to reactions
        mapped_reacs = map_ec_to_reac(mapped_reacs[['ec-code','ncbiprotein']])
        
        # @TODO the stuff below also appear multiple times
        # save for manual curation
        self.manual_curation['reacs, no mapping'] = mapped_reacs[mapped_reacs['id'].isnull()]
        self._statistics['reactions']['no NCBI, no EC'] = len(self.manual_curation['reacs, no mapping'])
        # map to model
        mapped_reacs = mapped_reacs[~mapped_reacs['id'].isnull()]
        mapped_reacs['add_to_GPR'] = mapped_reacs.apply(lambda x: self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1)
        
        # update attributes
        self.missing_genes = updated_missing_genes
        self.missing_reacs = mapped_reacs
    
    
    
    
############################################################################
# functions
############################################################################




############################################################################
# For filtering
############################################################################
# Evtl hier:
# Inspired by Dr. Reihaneh Mostolizadeh's function to add BioCyc reactions to a model
def replace_reaction_direction_with_fluxes(missing_reacs: pd.DataFrame) -> pd.DataFrame:
   """Extracts the flux lower & upper bounds for each reaction through the entries in column 'Reaction-Direction'
   
   Args:
      - missing_reacs (pd.DataFrame): 
         Table containing reactions & the respective Reaction-Directions
         
   Returns:
      pd.DataFrame: 
         Input table extended with the fluxes lower & upper bounds obtained from 
         the Reaction-Directions
   """
    
   def get_fluxes(row: pd.Series) -> dict[str: str]:
      direction = row['Reaction-Direction']
      fluxes = {}
      
      if type(direction) == float:
         # Use default bounds as described in readthedocs from COBRApy
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'RIGHT-TO-LEFT' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_0_bound'
      elif 'LEFT-TO-RIGHT' in direction:
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'REVERSIBLE' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_default_ub'
      else:
         #@TODO LOGGING.WARNING + Set to reversible
         pass
      
      return str(fluxes)
   
   missing_reacs['fluxes'] = missing_reacs.apply(get_fluxes, axis=1)
   missing_reacs.drop('Reaction-Direction', axis=1, inplace=True)
   
   return missing_reacs