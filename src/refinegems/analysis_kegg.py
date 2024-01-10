#!/usr/bin/env python
""" Provides functions to compare genes found in KEGG and in the model

Extracts all KEGG ids from the annotations and compares them to a list for you organism from KEGG.
Reactions with KEGG ids not found in the model are expanded to a table containing the KEGG id,
the locus_tag (old format), the EC number, the BiGG id and the Gene Protein Reaction (GPR) rule.
This section needs a gff file of your organism, the txt file from BiGG containing all reactions
and the KEGG identifier of your organism.

Due to the KEGG REST API this is relatively slow (model of size 1500 reactions - 20 min).
"""

import pandas as pd
from bioservices.kegg import KEGG
import io
import re
from libsbml import Model as libModel
from refinegems.io import parse_gff_for_gp_info
from refinegems.entities import get_model_genes, compare_gene_lists, get_model_reacs_or_metabs
from refinegems.analysis_db import get_bigg2other_db, compare_bigg_model

__author__ = "Famke Baeuerle"


def get_kegg_genes(organismid: str) -> pd.DataFrame:
    """Extracts list of genes from KEGG given an organism

    Args:
        - organismid (str): KEGG ID of organism which the model is based on

    Returns:
        pd.DataFrame: Table of all genes denoted in KEGG for the organism
    """
    k = KEGG()
    gene_list = k.list(organismid)

    return pd.read_table(io.StringIO(gene_list), header=None)


def get_locus_ec(genes_kegg_notmodel: pd.DataFrame) -> pd.DataFrame:
    """Creates columns with EC numbers for the locus tags of the genes

    Args:
        - genes_kegg_notmodel (pd.DataFrame): Genes present in KEGG but not in the model

    Returns:
        pd.DataFrame: Table of genes with locus tag and EC number
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


def get_locus_ec_kegg(locus_ec: pd.DataFrame) -> pd.DataFrame:
    """Searches for KEGG reactions based on EC numbers

    Args:
        - locus_ec (pd.DataFrame): Genes with locus tag and EC number

    Returns:
        pd.DataFrame: Table of genes with locus tag, EC number and KEGG Id
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


def get_locus_ec_kegg_bigg(locus_ec_kegg: pd.DataFrame, bigg_kegg: pd.DataFrame) -> pd.DataFrame:
    """Merges table with genes from model with BiGG / KEGG mapping to add BiGG Ids

    Args:
        - locus_ec_kegg (pd.DataFrame): Genes with locus tag, EC number and KEGG Id
        - bigg_kegg (pd.DataFrame): BiGG IDs with corresponding KEGG Ids

    Returns:
        pd.DataFrame: Table of genes with locus tag, EC number, KEGG Id and BiGG Id
    """
    locus_ec_kegg_bigg = locus_ec_kegg.merge(bigg_kegg, on=['KEGG'])
    return locus_ec_kegg_bigg


def get_locus_ec_kegg_bigg_gpr(locus_ec_kegg_bigg: pd.DataFrame, locus_gpr: pd.DataFrame) -> pd.DataFrame:
    """Merges table with genes from model if locus tag / GPR mapping to add GPRs

    Args:
        - locus_ec_kegg_bigg (pd.DataFrame): Genes with locus tag, EC number, KEGG Id and BiGG Id
        - locus_gpr (pd.DataFrame): Mapping from locus tags to GPRs

    Returns:
        pd.DataFrame: Table of genes with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """

    def slice_locus(locus):
        return locus[:-1]

    locus_ec_kegg_bigg['locus_tag'] = locus_ec_kegg_bigg.apply(
        lambda row: slice_locus(row['locus_tag']), axis=1)

    return locus_ec_kegg_bigg.merge(locus_gpr, how='left', on='locus_tag')


def kegg_gene_comp(model: libModel, organismid: str, gff_file: str) -> pd.DataFrame:
    """Exectues all steps to compare genes of the model to KEGG genes

    Args:
        - model (libModel): Model loaded with libSBML
        - organismid (str): KEGG ID of organism which the model is based on
        - gff_file (str): Path to gff file of organism of interest

    Returns:
        pd.DataFrame: Table containing missing reactions with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """
    model_genes = get_model_genes(model, True)
    model_reactions = get_model_reacs_or_metabs(model)
    kegg_genes = get_kegg_genes(organismid)
    bigg_kegg = get_bigg2other_db('KEGG')
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
