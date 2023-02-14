#!/usr/bin/env python
import pandas as pd
from Bio import Entrez
from libsbml import *
from refinegems.cvterms import add_cv_term_genes, add_cv_term_metabolites, add_cv_term_reactions
from refinegems.sboann import *
from refinegems.io import search_ncbi_for_gpr

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_model_genes(model: Model, kegg: bool=False) -> pd.DataFrame:
    """Extracts KEGG Genes/Locus tags from given model

    Args:
        model (model-libsbml): model loaded with libsbml
        kegg (bool): True if KEGG Genes should be extracted, otherwise False

    Returns:
        df: table with all KEGG Genes/Locus tags in the model
    """
    genes_in_model = []
    for gene in model.getPlugin(0).getListOfGeneProducts():
        if kegg:
            cv_terms = gene.getCVTerms()
            if cv_terms:
                for cv_term in cv_terms:
                    for idx in range(cv_term.getNumResources()):
                        uri = cv_term.getResourceURI(idx)
                        if 'kegg.genes' in uri: genes_in_model.append(uri.split('kegg.genes:')[1])
        else:
            genes_in_model.append(gene.getLabel())

    return pd.DataFrame(genes_in_model)


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_gene_lists(gps_in_model: pd.DataFrame, db_genes: pd.DataFrame, kegg: bool=True) -> pd.DataFrame:
    """Compares the provided dataframes according to column 0/'Locus_tag'
    
        Args:
            gps_in_model (DataFrame): pandas dataframe containing the KEGG Gene IDs/Locus tags in the model
            db_genes (DataFrame): pandas dataframe containing the KEGG Gene IDs for the organism from KEGG/
                                    locus tags (Accession-2) from BioCyc
    """
    in_db = db_genes.set_index(0) if kegg else db_genes.set_index('locus_tag')
    in_model = gps_in_model.set_index(0)
        
    genes_in_db_not_in_model = in_db[~in_db.index.isin(in_model.index)]
    
    return genes_in_db_not_in_model.reset_index().iloc[:, 0] if kegg else genes_in_db_not_in_model.reset_index()


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_model_reacs_or_metabs(model_libsbml: Model, metabolites: bool=False):
    """Extracts table of reactions/metabolites with BiGG Ids from model

    Args:
        model_libsbml (Model): model loaded with libsbml
        metabolites (bool): set to True if metabolites from model should be extracted

    Returns:
        df: table with BiGG Ids of reactions in the model
    """
    reac_or_metab_list = model_libsbml.getListOfSpecies() if metabolites else model_libsbml.getListOfReactions()

    list_of_reacs_or_metabs = []
    for reac_or_metab in reac_or_metab_list:
        list_of_reacs_or_metabs.append(reac_or_metab.id[2:])

    reac_or_metab_list_df = pd.Series(list_of_reacs_or_metabs)
    reac_or_metab_list_df = pd.DataFrame(reac_or_metab_list_df, columns=['bigg_id'])

    return reac_or_metab_list_df


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
    name, locus = search_ncbi_for_gpr(locus_tag)
    gpr = model.getPlugin(0).createGeneProduct()
    gpr.setName(name)
    gpr.setId(locus_tag)
    gpr.setMetaId('meta_' + locus_tag)
    gpr.setLabel(locus_tag)
    gpr.setSBOTerm("SBO:0000243")
    add_cv_term_genes(locus_tag, 'NCBI', gpr)
    return gpr, model


def create_gp(model: Model, model_id: str, name: str, locus_tag: str, protein_id: str) -> tuple[GeneProduct, Model]:
    """creates GeneProduct in the given model

    Args:
        model (libsbml-model): model loaded with libSBML
        model_id (str): ID identical to ID that CarveMe adds from the NCBI FASTA input file
        name (str): Name of the GeneProduct
        locus_tag (str): genome-specific locus tag used as label in the model
        protein_id (str): NCBI Protein/RefSeq ID

    Returns:
        tuple: (gpr, modified model)
    """
    gp = model.getPlugin(0).createGeneProduct()
    gp.setId(model_id)
    gp.setName(name)
    gp.setLabel(locus_tag)
    gp.setSBOTerm('SBO:0000243')
    gp.setMetaId(f'meta_{model_id}')
    add_cv_term_genes(protein_id, 'NCBI', gp)
    return gp, model


def create_species(model: Model, metabolite_id: str, name: str, compartment_id: str, charge: int, chem_formula: str):
    """creates Species/Metabolite in the given model

    Args:
        model (libsbml-model): model loaded with libSBML
        metabolite_id (string): metabolite ID within model (If model from CarveMe, preferable a BiGG ID)
        name (str): name of the metabolite
        compartment_id (str): Id of the compartment where metabolite resides
        chem_formula (str): chemical formula for the metabolite

    Returns:
        tuple: (metabolite, modified model)
    """
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
    metabolite.getPlugin(0).setCharge(charge)
    metabolite.getPlugin(0).setChemicalFormula(chem_formula)
    add_cv_term_metabolites(metabolite_id, 'BIGG', metabolite)
    return metabolite, model


def get_fluxes() -> dict[str: int]:
    # return dict of fluxes
    fluxes = {}
    return fluxes


def get_reversible():
    pass


def create_reaction(model, reaction_id, name, reactants, products, fluxes, reversible, fast, sbo=None):
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
    sbo = sbo if sbo else 'SBO:0000167'  # SBO term for biochemical or transport reaction
    reaction.setSBOTerm(sbo)
    fast = fast if fast else False
    reaction.setFast(fast)
    reversible = reversible if reversible else get_reversible()
    reaction.setReversible(reversible)
    for metab, stoich in reactants.items(): #reactants as dict with metab:stoich
        reaction.addReactant(model.getSpecies('M_' + metab), stoich)
    for metab, stoich in products.items(): #reactants as dict with metab:stoich
        reaction.addProduct(model.getSpecies('M_' + metab), stoich)
    fluxes = fluxes if fluxes else get_fluxes()
    reaction.getPlugin(0).setLowerFluxBound(fluxes['lower_bound'])
    reaction.getPlugin(0).setUpperFluxBound(fluxes['upper_bound'])
    add_cv_term_reactions(reaction_id, 'BIGG', reaction)
    return reaction, model
 