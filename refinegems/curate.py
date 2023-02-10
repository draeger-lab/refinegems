#!/usr/bin/env python
""" Functions to enable annotation of entities using a manual curated table

While working on GEMs the user might come across ill-annotated or missing metabolites, reactions and genes. This module aims to enable faster manual curation by  allowing to edit an excel table directly which is used to update the given model. This module makes use of the cvterms module aswell.
"""
from tqdm.auto import tqdm
from refinegems.cvterms import add_cv_term_reactions, add_cv_term_metabolites, metabol_db_dict, parse_id_from_cv_term
from refinegems.entities import create_gpr, create_reaction

__author__ = "Famke Baeuerle"
    
            
def add_reactions_from_table(model, table, email):
    """adds all reactions with their info given in the table to the given model
       (wrapper function to use with table given in the data folder)

    Args:
        model (libsbml-model): model loaded with libSBML
        table (pd-DataFrame): manual_curation tabel loaded with pandas, sheet gapfill
        email (string): User Email to access the NCBI Entrez database

    Returns:
        model: modified model with new reactions
    """
    for reaction_info in tqdm(table.groupby('BIGG')):
        reac_id = reaction_info[0]
        if model.getReaction(str(reac_id)) is None:
            reactants = dict(table.loc[table['BIGG'] == reac_id, ['educts', 'stoich_e']].dropna().values)
            products = dict(table.loc[table['BIGG'] == reac_id, ['products', 'stoich_p']].dropna().values)
            fluxes = table.loc[table['BIGG'] == reac_id, ['lower_bound', 'upper_bound']].dropna().to_dict('records')[0]
            name = table.loc[table['BIGG'] == reac_id, ['name']].dropna().to_dict('records')[0]['name']
            reversible = table.loc[table['BIGG'] == reac_id, ['reversible']].dropna().to_dict('records')[0]['reversible']
            fast = table.loc[table['BIGG'] == reac_id, ['fast']].dropna().to_dict('records')[0]['fast']
            try:
                sbo = table.loc[table['BIGG'] == reac_id, ['sbo']].dropna().to_dict('records')[0]['sbo']
            except (IndexError):
                print('SBO Term for ' + str(reac_id) + ' will be set to standard "SBO:0000167" (biochemical or transport reaction)')
                sbo = "SBO:0000167"
            reaction, model = create_reaction(model, reac_id, name, reactants, products, fluxes, reversible, fast, sbo)
            for (columnName, columnData) in table.loc[table['BIGG'] == reac_id].drop(['educts', 'stoich_e','products', 'stoich_p','lower_bound', 'upper_bound', 'name', 'reversible', 'fast'], axis=1).fillna(0).iteritems():
                for entry in columnData.values:
                    if not entry == 0:
                        if columnName == 'locus':
                            reaction.getPlugin(0).createGeneProductAssociation().createGeneProductRef().setGeneProduct(str(entry))
                            if model.getPlugin(0).getGeneProductByLabel(str(entry)) is None:
                                gpr, model = create_gpr(model, str(entry), email)
                        else:
                            add_cv_term_reactions(str(entry), str(columnName), reaction)
    return model


def update_annotations_from_table(model, table):
    """updates annotation of metabolites given in the table
       (wrapper function to use with table given in the data folder)

    Args:
        model (libsbml-model): model loaded with libSBML
        table (pd-DataFrame): manual_curation tabel loaded with pandas, sheet metabs

    Returns:
        model: modified model with new annotations
    """
    table = table.drop(['Name', 'FORMULA', 'Notiz'], axis=1).fillna(0)
    table['PUBCHEM'] = table['PUBCHEM'].astype(int)
    for metab_info in tqdm(table.groupby('BIGG')):
        met = metab_info[0]
        for comp in ['_c', '_e', '_p']:
            try:
                metab = model.getSpecies('M_' + met + comp)
                #metab.unsetAnnotation()
                if not metab.isSetMetaId():
                    metab.setMetaId('meta_' + metab.getId())
                for (columnName, columnData) in table.loc[table['BIGG'] == met].iteritems():
                    for entry in columnData.values:
                        if not entry == 0:
                            add_cv_term_metabolites(str(entry), str(columnName), metab)
            except (AttributeError):
                print(met + comp + ' not in model')
    return model

def update_annotations_from_others(model):
    """Synchronizes metabolite annotations for core, periplasm and extracelullar

    Args:
        model (libsbml-model): model loaded with libSBML

    Returns:
        model: modified with synchronized annotations
    """
    for metab in model.getListOfSpecies():
        base = metab.getId()[:-2]
        for comp in ['_c', '_e', '_p']:
            other_metab = model.getSpecies(base + comp)
            if other_metab is not None:
                if not other_metab.isSetMetaId():
                    other_metab.setMetaId('meta_' + other_metab.getId())
                for db_id, code in metabol_db_dict.items():
                    id = parse_id_from_cv_term(metab, code)
                    for entry in id:
                        if entry is not None:
                            add_cv_term_metabolites(entry, db_id, other_metab)    
    return model               