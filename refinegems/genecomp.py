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
import io, re
import gffutils

__author__ = "Famke Baeuerle"

def get_model_genes(model):
    genes_in_model = []
    for gene in model.getPlugin(0).getListOfGeneProducts():
        cv_terms = gene.getCVTerms()
        if cv_terms:
            for cv_term in cv_terms:
                for idx in range(cv_term.getNumResources()):
                    uri = cv_term.getResourceURI(idx)
                    genes_in_model.append(uri[-16:])
                    
    return pd.DataFrame(genes_in_model)

def get_kegg_genes(organismid):
    k = KEGG()
    gene_list = k.list(organismid)

    return pd.read_table(io.StringIO(gene_list), header=None)

def compare_kegg_model(model_genes, kegg_genes):
    ke = kegg_genes.set_index(0)
    mo = model_genes.set_index(0)
    genes_in_kegg_not_in_model = ke[~ke.index.isin(mo.index)]
    return genes_in_kegg_not_in_model.reset_index().iloc[:,0]

def get_locus_ec(genes_kegg_notmodel):
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
                real_ec[entry[:12]] = '['+m.group(0)
                
    locus_ec = pd.DataFrame.from_dict(real_ec, orient='index').reset_index().rename(columns={'index':'locus_tag', 0:'EC-number'})
    
    def slice_ec(ec):
        new = ec[4:]
        new2 = new[:-1]
        return new2

    locus_ec['EC'] = locus_ec.apply(lambda row: slice_ec(row['EC-number']), axis=1)
    locus_ec = locus_ec.drop('EC-number', axis=1) 
    
    return locus_ec

def get_locus_ec_kegg(locus_ec):
    
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
    
    locus_ec['KEGG_Ids'] = locus_ec.apply(lambda row: get_kegg_reaction(row['EC']), axis=1)
    locus_ec = locus_ec.dropna()
    locus_ec['KEGG_Ids2'] = locus_ec.apply(lambda row: drop_nonreac(row['KEGG_Ids']), axis=1)
    locus_ec = locus_ec.dropna()
    locus_ec['KEGG'] = locus_ec.apply(lambda row: slice_kegg(row['KEGG_Ids2']), axis=1)
    locus_ec_kegg = locus_ec.dropna().drop('KEGG_Ids', axis=1).drop('KEGG_Ids2', axis=1)
    
    return locus_ec_kegg

def get_bigg_kegg(biggreactions):
    # make the download of biggreactions possible to maintain database
    all_reac_bigg = pd.read_csv(biggreactions, sep='\t').drop('model_list', axis=1).dropna()
    
    def find_kegg(database_links):
        m = re.search('KEGG Reaction: http://identifiers.org/kegg.reaction/(.*?); ', database_links)
        if m:
            return m.group(0)
        else:
            return None
        
    def slice_kegg(kegg):
        return kegg[52:-2]
    
    all_reac_bigg['KEGG_link'] = all_reac_bigg.apply(lambda row: find_kegg(row['database_links']), axis=1)
    all_reac_bigg['KEGG'] = all_reac_bigg.dropna().apply(lambda row: slice_kegg(row['KEGG_link']), axis=1)
    all_reac_bigg = all_reac_bigg.dropna().drop('KEGG_link', axis=1)

    return all_reac_bigg[['bigg_id','KEGG']]

def get_locus_gpr(gff_file):
    db = gffutils.create_db(gff_file, ':memory:',  merge_strategy='create_unique')
    mapping_cds = {}
    for feature in db.all_features():
        attr = dict(feature.attributes)
        try: 
            if str(attr['gbkey'][0]) == 'CDS':
                mapping_cds[attr['Name'][0]] = attr['Parent'][0]
        except:
            pass
    mapping_df = pd.DataFrame.from_dict(mapping_cds, columns=['Parent'], orient='index').reset_index().rename(columns={'index':'GPR'})

    def extract_locus(feature):
        try:
            return db[feature].attributes['old_locus_tag'][0]
        except:
            pass
        return None

    mapping_df['locus_tag'] = mapping_df.apply(lambda row: extract_locus(row['Parent']), axis=1)
    return mapping_df.drop('Parent', axis=1)

def get_locus_ec_kegg_bigg(locus_ec_kegg, bigg_kegg):
    locus_ec_kegg_bigg = locus_ec_kegg.merge(bigg_kegg, on=['KEGG'])
    return locus_ec_kegg_bigg

def get_locus_ec_kegg_bigg_gpr(locus_ec_kegg_bigg, locus_gpr):

    def slice_locus(locus):
        return locus[:-1]

    locus_ec_kegg_bigg['locus_tag']=locus_ec_kegg_bigg.apply(lambda row: slice_locus(row['locus_tag']), axis=1)

    return locus_ec_kegg_bigg.merge(locus_gpr, how='left', on='locus_tag')

def get_model_reactions(model):
    reac_list = model.getListOfReactions()

    list_of_reac = []
    for reac in reac_list:
        list_of_reac.append(reac.id[2:])

    reac_list_df = pd.Series(list_of_reac)
    reac_list_df = pd.DataFrame(reac_list_df, columns=['bigg_id'])

    return reac_list_df

def compare_bigg_model(locus_ec_kegg_bigg_gpr, model_reactions):
    mapp = locus_ec_kegg_bigg_gpr.set_index('bigg_id')
    reacs = model_reactions.set_index('bigg_id')

    reactions_missing_in_model = mapp[~mapp.index.isin(reacs.index)].reset_index()

    ambig_kegg = locus_ec_kegg_bigg_gpr.loc[locus_ec_kegg_bigg_gpr.duplicated(subset=['KEGG'], keep=False)]
    ambig_kegg = ambig_kegg.set_index('KEGG').drop(['locus_tag', 'EC'], axis=1).sort_index()

    ambig = ambig_kegg.set_index('bigg_id')
    miss = reactions_missing_in_model.set_index('bigg_id')

    reactions_missing_in_model_non_dup = miss[~miss.index.isin(ambig.index)].reset_index()

    return reactions_missing_in_model_non_dup

def genecomp(model, organismid, biggreactions, gff_file):
    model_genes = get_model_genes(model)
    model_reactions = get_model_reactions(model)
    kegg_genes = get_kegg_genes(organismid)
    bigg_kegg = get_bigg_kegg(biggreactions)
    genes_kegg_notmodel = compare_kegg_model(model_genes, kegg_genes)
    locus_gpr = get_locus_gpr(gff_file)
    locus_ec = get_locus_ec(genes_kegg_notmodel)
    locus_ec_kegg = get_locus_ec_kegg(locus_ec)
    locus_ec_kegg_bigg = get_locus_ec_kegg_bigg(locus_ec_kegg, bigg_kegg)
    locus_ec_kegg_bigg_gpr = get_locus_ec_kegg_bigg_gpr(locus_ec_kegg_bigg, locus_gpr)
    missing_reactions = compare_bigg_model(locus_ec_kegg_bigg_gpr, model_reactions)
    return missing_reactions