#!/usr/bin/env python
""" Functions to enable annotation of entities using a manual curated table

While working on GEMs the user might come across ill-annotated or missing metabolites, reactions and genes. This module aims to enable faster manual curation by  allowing to edit an excel table directly which is used to update the given model. This module makes use of the cvterms module aswell.
"""
from Bio import Entrez, SeqIO
from tqdm.auto import tqdm
from refinegems.cvterms import add_cv_term_reactions, add_cv_term_genes, add_cv_term_metabolites

__author__ = "Famke Baeuerle"


def get_name_from_locus(locus):
    """fetches protein name from NCBI

    Args:
        locus (string): NCBI compatible locus_tag

    Returns:
        str: protein name / description
    """
    handle = Entrez.efetch(
        db="protein",
        id=locus,
        rettype="gbwithparts",
        retmode='text')
    records = SeqIO.parse(handle, "gb")

    for i, record in enumerate(records):
        return record.description

def create_reaction(model, reaction_id, name, reactants, products, fluxes, reversible, fast, sbo):
    """creates new reaction in the given model

    Args:
        model (libsbml-model): model loaded with libSBML
        reaction_id (str): BiGG id of the reaction to create
        name (string): human readable name of the reaction
        reactants (dict): metabolites as keys and their stoichiometry as values
        products (dict): metabolites as keys and their stoichiometry as values
        fluxes (dict): lower_bound and upper_bound as keys
        reversible (bool): true/false for the reaction
        fast (bool): true/false for the reaction
        sbo (string): SBO term of the reaction

    Returns:
        tuple: (reaction, modified model)
    """
    reaction = model.createReaction()
    reaction.setId('R_' + reaction_id)
    reaction.setName(name)
    reaction.setMetaId('meta_R_' + reaction_id)
    reaction.setSBOTerm(sbo)
    reaction.setFast(fast)
    reaction.setReversible(reversible)
    for metab, stoich in reactants.items(): #reactants as dict with metab:stoich
        reaction.addReactant(model.getSpecies('M_' + metab), stoich)
    for metab, stoich in products.items(): #reactants as dict with metab:stoich
        reaction.addProduct(model.getSpecies('M_' + metab), stoich)
    reaction.getPlugin(0).setLowerFluxBound(fluxes['lower_bound'])
    reaction.getPlugin(0).setUpperFluxBound(fluxes['upper_bound'])
    return reaction, model
    
def create_gpr(model, locus_tag, email):
    """creates GeneProduct in the given model

    Args:
        model (libsbml-model): model loaded with libSBML
        locus_tag (string): NCBI compatible locus_tag
        email (string): User Email to access the NCBI Entrez database

    Returns:
        tuple: (gpr, modified model)
    """
    Entrez.email = email
    gpr = model.getPlugin(0).createGeneProduct()
    name = get_name_from_locus(locus_tag)
    gpr.setName(name)
    gpr.setId(locus_tag)
    gpr.setMetaId('meta_' + locus_tag)
    gpr.setLabel(locus_tag)
    gpr.setSBOTerm("SBO:0000243")
    add_cv_term_genes(locus_tag, 'NCBI', gpr)
    return gpr, model
            
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
                metab.unsetAnnotation()
                for (columnName, columnData) in table.loc[table['BIGG'] == met].iteritems():
                    for entry in columnData.values:
                        if not entry == 0:
                            add_cv_term_metabolites(str(entry), str(columnName), metab)
            except (AttributeError):
                print(met + comp + ' not in model')
    return model