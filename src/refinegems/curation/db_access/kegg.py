#!/usr/bin/env python
"""Retrieve and parse information from the KEGG database. 

old text (most functionality described here has been moved or deprecated)

Provides functions to compare genes found in KEGG and in the model

Extracts all KEGG ids from the annotations and compares them to a list for you organism from KEGG.
Reactions with KEGG ids not found in the model are expanded to a table containing the KEGG id,
the locus_tag (old format), the EC number, the BiGG id and the Gene Protein Reaction (GPR) rule.
This section needs a gff file of your organism, the txt file from BiGG containing all reactions
and the KEGG identifier of your organism.

Due to the KEGG REST API this is relatively slow (model of size 1500 reactions - 20 min).
"""

__author__ = "Famke Baeuerle and Carolin Brune"

############################################################################
# requirements
############################################################################

import io
import pandas as pd
import re

from bioservices.kegg import KEGG
from Bio.KEGG import REST, Gene, Enzyme
from libsbml import Model as libModel

from ...utility.io import parse_gff_for_gp_info
from ...utility.entities import get_model_genes, compare_gene_lists, get_model_reacs_or_metabs
from .db import get_bigg_db_mapping, compare_bigg_model

############################################################################
# functions
############################################################################


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
            The locus in the format <orgnismid>:<locus_tag>

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
        # gett KEGG orthology ID
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


# re-check for deprecation
# ------------------------

# @DEPRECATED : still needed ?
def get_locus_ec(genes_kegg_notmodel: pd.DataFrame) -> pd.DataFrame:
    """Creates columns with EC numbers for the locus tags of the genes

    Args:
        - genes_kegg_notmodel (pd.DataFrame): 
            Genes present in KEGG but not in the model

    Returns:
        pd.DataFrame: 
            Table of genes with locus tag and EC number
    """
    k = KEGG()

    ec_dict = {}
    for gene in genes_kegg_notmodel:
        entry = k.parse(k.get(gene))
        try:
            ec_dict[entry['ENTRY']] = (entry['ORTHOLOGY'])
        except(KeyError):
            pass

    real_ec = {}
    for entry, ortho in ec_dict.items():
        for key, value in ortho.items():
            m = re.search('(?:EC).*', value)
            if m:
                real_ec[entry[:12]] = '[' + m.group(0)

    locus_ec = pd.DataFrame.from_dict(
        real_ec,
        orient='index').reset_index().rename(
        columns={
            'index': 'locus_tag',
            0: 'EC-number'})

    def slice_ec(ec):
        new = ec[4:]
        new2 = new[:-1]
        return new2

    locus_ec['EC'] = locus_ec.apply(
        lambda row: slice_ec(
            row['EC-number']), axis=1)
    locus_ec = locus_ec.drop('EC-number', axis=1)

    return locus_ec


# @DEPRECATED : still needed ?
def get_locus_ec_kegg(locus_ec: pd.DataFrame) -> pd.DataFrame:
    """Searches for KEGG reactions based on EC numbers

    Args:
        - locus_ec (pd.DataFrame): 
            Genes with locus tag and EC number

    Returns:
        pd.DataFrame: 
            Table of genes with locus tag, EC number and KEGG Id
    """

    def get_kegg_reaction(ec_number):
        k = KEGG()
        gene = k.parse(k.get(ec_number))
        try:
            return gene['REACTION'][-1]
        except(KeyError):
            pass
        return None

    def drop_nonreac(kegg_id):
        if len(kegg_id) != 11:
            return None
        else:
            return kegg_id

    def slice_kegg(kegg):
        return kegg[4:-1]

    locus_ec['KEGG_Ids'] = locus_ec.apply(
        lambda row: get_kegg_reaction(
            row['EC']), axis=1)
    locus_ec = locus_ec.dropna()
    locus_ec.loc[:, 'KEGG_Ids2'] = locus_ec.apply(
        lambda row: drop_nonreac(
            row['KEGG_Ids']), axis=1)
    locus_ec = locus_ec.dropna()
    locus_ec['KEGG'] = locus_ec.apply(
        lambda row: slice_kegg(
            row['KEGG_Ids2']), axis=1)
    locus_ec_kegg = locus_ec.dropna().drop(
        'KEGG_Ids', axis=1).drop(
        'KEGG_Ids2', axis=1)

    return locus_ec_kegg


# @DEPRECATED : still needed ?
def get_locus_ec_kegg_bigg(locus_ec_kegg: pd.DataFrame, bigg_kegg: pd.DataFrame) -> pd.DataFrame:
    """Merges table with genes from model with BiGG / KEGG mapping to add BiGG Ids

    Args:
        - locus_ec_kegg (pd.DataFrame): 
            Genes with locus tag, EC number and KEGG Id
        - bigg_kegg (pd.DataFrame): 
            BiGG IDs with corresponding KEGG Ids

    Returns:
        pd.DataFrame: 
            Table of genes with locus tag, EC number, KEGG Id and BiGG Id
    """
    locus_ec_kegg_bigg = locus_ec_kegg.merge(bigg_kegg, on=['KEGG'])
    return locus_ec_kegg_bigg


# @DEPRECATED : still needed ?
def get_locus_ec_kegg_bigg_gpr(locus_ec_kegg_bigg: pd.DataFrame, locus_gpr: pd.DataFrame) -> pd.DataFrame:
    """Merges table with genes from model if locus tag / GPR mapping to add GPRs

    Args:
        - locus_ec_kegg_bigg (pd.DataFrame): 
            Genes with locus tag, EC number, KEGG Id and BiGG Id
        - locus_gpr (pd.DataFrame): 
            Mapping from locus tags to GPRs

    Returns:
        pd.DataFrame: 
            Table of genes with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """

    def slice_locus(locus):
        return locus[:-1]

    locus_ec_kegg_bigg['locus_tag'] = locus_ec_kegg_bigg.apply(
        lambda row: slice_locus(row['locus_tag']), axis=1)

    return locus_ec_kegg_bigg.merge(locus_gpr, how='left', on='locus_tag')


# @DEPRECATED - basically an old version of first part of the KEGGapFiller
def kegg_gene_comp(model: libModel, organismid: str, gff_file: str) -> pd.DataFrame:
    """Exectues all steps to compare genes of the model to KEGG genes

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - organismid (str):
            KEGG ID of organism which the model is based on
        - gff_file (str): 
            Path to gff file of organism of interest

    Returns:
        pd.DataFrame: 
            Table containing missing reactions with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """
    model_genes = get_model_genes(model, True)
    model_reactions = get_model_reacs_or_metabs(model)
    kegg_genes = get_kegg_genes(organismid)
    bigg_kegg = get_bigg_db_mapping('KEGG',False)
    genes_kegg_notmodel = compare_gene_lists(model_genes, kegg_genes)
    locus_gpr = parse_gff_for_gp_info(gff_file)
    locus_ec = get_locus_ec(genes_kegg_notmodel)
    locus_ec_kegg = get_locus_ec_kegg(locus_ec)
    locus_ec_kegg_bigg = get_locus_ec_kegg_bigg(locus_ec_kegg, bigg_kegg)
    locus_ec_kegg_bigg_gpr = get_locus_ec_kegg_bigg_gpr(
        locus_ec_kegg_bigg, locus_gpr)
    missing_reactions = compare_bigg_model(
        locus_ec_kegg_bigg_gpr, model_reactions)
    return missing_reactions

