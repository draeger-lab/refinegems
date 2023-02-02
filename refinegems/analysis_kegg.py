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
import gffutils
from refinegems.gapfill import get_model_genes, compare_gene_lists

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


def get_kegg_genes(organismid):
    """Extracts list of genes from KEGG given an organism

    Args:
        organismid (Str): KEGG Id of organism which the model is based on

    Returns:
        df: table of all genes denoted in KEGG for the organism
    """
    k = KEGG()
    gene_list = k.list(organismid)

    return pd.read_table(io.StringIO(gene_list), header=None)


def get_locus_ec(genes_kegg_notmodel):
    """Creates columns with EC numbers for the locus tags of the genes

    Args:
        genes_kegg_notmodel (df): genes present in KEGG but not in the model

    Returns:
        df: table of genes with locus tag and EC number
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


def get_locus_ec_kegg(locus_ec):
    """Searches for KEGG reactions based on EC numbers

    Args:
        locus_ec (df): genes with locus tag and EC number

    Returns:
        df: table of genes with locus tag, EC number and KEGG Id
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
    locus_ec['KEGG_Ids2'] = locus_ec.apply(
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


def get_bigg_kegg(biggreactions):
    """Uses list of BiGG reactions to get a mapping from BiGG to KEGG Id

    Args:
        biggreactions (Str): path to file containing BiGG database

    Returns:
        df: table containing BiGG Ids with corresponding KEGG Ids
    """
    # make the download of biggreactions possible to maintain database
    all_reac_bigg = pd.read_csv(
        biggreactions,
        sep='\t').drop(
        'model_list',
        axis=1).dropna()

    def find_kegg(database_links):
        m = re.search(
            'KEGG Reaction: http://identifiers.org/kegg.reaction/(.*?); ',
            database_links)
        if m:
            return m.group(0)
        else:
            return None

    def slice_kegg(kegg):
        return kegg[52:-2]

    all_reac_bigg['KEGG_link'] = all_reac_bigg.apply(
        lambda row: find_kegg(row['database_links']), axis=1)
    all_reac_bigg['KEGG'] = all_reac_bigg.dropna().apply(
        lambda row: slice_kegg(row['KEGG_link']), axis=1)
    all_reac_bigg = all_reac_bigg.dropna().drop('KEGG_link', axis=1)

    return all_reac_bigg[['bigg_id', 'KEGG']]


def get_locus_gpr(gff_file):
    """Searches gff file of organism for gene protein reactions based on locus tags

    Args:
        gff_file (Str): path to gff file of organism of interest

    Returns:
        df: table containing mapping from locus tag to GPR
    """
    db = gffutils.create_db(
        gff_file,
        ':memory:',
        merge_strategy='create_unique')
    mapping_cds = {}
    for feature in db.all_features():
        attr = dict(feature.attributes)
        try:
            if str(attr['gbkey'][0]) == 'CDS':
                mapping_cds[attr['Name'][0]] = attr['Parent'][0]
        except BaseException:
            pass
    mapping_df = pd.DataFrame.from_dict(
        mapping_cds,
        columns=['Parent'],
        orient='index').reset_index().rename(
        columns={
            'index': 'GPR'})

    def extract_locus(feature):
        try:
            return db[feature].attributes['old_locus_tag'][0]
        except BaseException:
            pass
        return None

    mapping_df['locus_tag'] = mapping_df.apply(
        lambda row: extract_locus(row['Parent']), axis=1)
    return mapping_df.drop('Parent', axis=1)


def get_locus_ec_kegg_bigg(locus_ec_kegg, bigg_kegg):
    """Merges table with genes from model with BiGG / KEGG mapping to add BiGG Ids

    Args:
        locus_ec_kegg (df): genes with locus tag, EC number and KEGG Id
        bigg_kegg (df): BiGG Ids with corresponding KEGG Ids

    Returns:
        df: table of genes with locus tag, EC number, KEGG Id and BiGG Id
    """
    locus_ec_kegg_bigg = locus_ec_kegg.merge(bigg_kegg, on=['KEGG'])
    return locus_ec_kegg_bigg


def get_locus_ec_kegg_bigg_gpr(locus_ec_kegg_bigg, locus_gpr):
    """Merges table with genes from model if locus tag / GPR mapping to add GPRs

    Args:
        locus_ec_kegg_bigg (df): genes with locus tag, EC number, KEGG Id and BiGG Id
        locus_gpr (df): mapping from locus tags to GPRs

    Returns:
        df: table of genes with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """

    def slice_locus(locus):
        return locus[:-1]

    locus_ec_kegg_bigg['locus_tag'] = locus_ec_kegg_bigg.apply(
        lambda row: slice_locus(row['locus_tag']), axis=1)

    return locus_ec_kegg_bigg.merge(locus_gpr, how='left', on='locus_tag')


def get_model_reactions(model):
    """Extracts table of reactions with BiGG Ids from model

    Args:
        model (libsbml-model): model loaded with libsbml

    Returns:
        df: table with BiGG Ids of reactions in the model
    """
    reac_list = model.getListOfReactions()

    list_of_reac = []
    for reac in reac_list:
        list_of_reac.append(reac.id[2:])

    reac_list_df = pd.Series(list_of_reac)
    reac_list_df = pd.DataFrame(reac_list_df, columns=['bigg_id'])

    return reac_list_df


def compare_bigg_model(locus_ec_kegg_bigg_gpr, model_reactions):
    """Compares reactions by genes extracted via KEGG to reactions in the model
        Needed to back check previous comparisons.

    Args:
        locus_ec_kegg_bigg_gpr (df): locus tag, EC number, KEGG Id, BiGG Id and GPR
        model_reactions (df): BiGG Ids of reactions in the model

    Returns:
        df: table containing reactions present in KEGG but not in the model
    """
    mapp = locus_ec_kegg_bigg_gpr.set_index('bigg_id')
    reacs = model_reactions.set_index('bigg_id')

    reactions_missing_in_model = mapp[~mapp.index.isin(
        reacs.index)].reset_index()

    ambig_kegg = locus_ec_kegg_bigg_gpr.loc[locus_ec_kegg_bigg_gpr.duplicated(
        subset=['KEGG'], keep=False)]
    ambig_kegg = ambig_kegg.set_index('KEGG').drop(
        ['locus_tag', 'EC'], axis=1).sort_index()

    ambig = ambig_kegg.set_index('bigg_id')
    miss = reactions_missing_in_model.set_index('bigg_id')

    reactions_missing_in_model_non_dup = miss[~miss.index.isin(
        ambig.index)].reset_index()

    return reactions_missing_in_model_non_dup


def kegg_gene_comp(model, organismid, biggreactions, gff_file):
    """Exectues all steps to compare genes of the model to KEGG genes

    Args:
        model (libsbml-model): model loaded with libsbml
        organismid (Str): KEGG Id of organism which the model is based on
        biggreactions (Str): path to file containing BiGG database
        gff_file (Str): path to gff file of organism of interest

    Returns:
        df: table containing missing reactions with locus tag, EC number, KEGG Id, BiGG Id and GPR
    """
    model_genes = get_model_genes(model, True)
    model_reactions = get_model_reactions(model)
    kegg_genes = get_kegg_genes(organismid)
    bigg_kegg = get_bigg_kegg(biggreactions)
    genes_kegg_notmodel = compare_gene_lists(model_genes, kegg_genes, True)
    locus_gpr = get_locus_gpr(gff_file)
    locus_ec = get_locus_ec(genes_kegg_notmodel)
    locus_ec_kegg = get_locus_ec_kegg(locus_ec)
    locus_ec_kegg_bigg = get_locus_ec_kegg_bigg(locus_ec_kegg, bigg_kegg)
    locus_ec_kegg_bigg_gpr = get_locus_ec_kegg_bigg_gpr(
        locus_ec_kegg_bigg, locus_gpr)
    missing_reactions = compare_bigg_model(
        locus_ec_kegg_bigg_gpr, model_reactions)
    return missing_reactions
