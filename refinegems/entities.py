#!/usr/bin/env python
from Bio import Entrez, SeqIO
from refinegems.cvterms import add_cv_term_genes, add_cv_term_metabolites, add_cv_term_reactions


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
    add_cv_term_reactions(reaction_id, 'BIGG', reaction)
    return reaction, model
 
 
def create_species(model: Model, metabolite_id: str, name: str, compartment_id: str):
   metabolite = model.createSpecies()
   metabolite.setId(f'M_{metabolite_id}')
   metabolite.setName(name)
   metabolite.setMetaId(f'meta_M_{metabolite_id}')
   metabolite.setSBOTerm('SBO:0000247')
   metabolite.setInitialAmount(float('NaN'))
   metabolite.setHasOnlySubstanceUnits(True)
   metabolite.setBoundaryCondition(False)
   metabolite.setConstant(False)
   metabolite.setCompartment(compartment_id)
   metabolite.setCharge(0)  # Reihaneh set charges to 0 as not provided by BioCyc
   metabolite.setChemicalFormula()
   metabolite.enablePackage()
   add_cv_term_metabolites(metabolite_id, 'BIGG', metabolite)
 