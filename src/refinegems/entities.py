#!/usr/bin/env python
import re
import pandas as pd
from Bio import Entrez
from libsbml import Model as libModel
from libsbml import GeneProduct, Species, Reaction
from refinegems.cvterms import add_cv_term_genes, add_cv_term_metabolites, add_cv_term_reactions
from refinegems.io import search_ncbi_for_gpr
from typing import Union

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel"


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_model_genes(model: libModel, kegg: bool=False) -> pd.DataFrame:
    """Extracts KEGG Genes/Locus tags from given model

    Args:
        - model (model-libsbml): Model loaded with libSBML
        - kegg (bool): True if KEGG Genes should be extracted, otherwise False

    Returns:
        pd.DataFrame: Table with all KEGG Genes/Locus tags in the model
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
    """Compares the provided tables according to column 0/'Locus_tag'
    
    Args:
        - gps_in_model (pd.DataFrame): Table containing the KEGG Gene IDs/Locus tags in the model
        - db_genes (pd.DataFrame): Table containing the KEGG Gene IDs for the organism from KEGG/
                                locus tags (Accession-2) from BioCyc
        - kegg (bool): True if KEGG Genes should be extracted, otherwise False
        
    Returns:
        pd.DataFrame: Table containing all missing genes
    """
    in_db = db_genes.set_index(0) if kegg else db_genes.set_index('locus_tag')
    in_model = gps_in_model.set_index(0)
        
    genes_in_db_not_in_model = in_db[~in_db.index.isin(in_model.index)]
    
    return genes_in_db_not_in_model.reset_index().iloc[:, 0] if kegg else genes_in_db_not_in_model.reset_index()


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_model_reacs_or_metabs(model_libsbml: libModel, metabolites: bool=False, col_name: str='bigg_id') -> pd.DataFrame:
    """Extracts table of reactions/metabolites with BiGG IDs from model

    Args:
        - model_libsbml (libModel): Model loaded with libSBML
        - metabolites (bool): Set to True if metabolites from model should be extracted
        - col_name (str): Name to be used for column in Table, default: 'bigg_id'

    Returns:
        pd.DataFrame: Table with model identifiers for either metabolites or reactions
    """
    reac_or_metab_list = model_libsbml.getListOfSpecies() if metabolites else model_libsbml.getListOfReactions()

    list_of_reacs_or_metabs = []
    for reac_or_metab in reac_or_metab_list:
        list_of_reacs_or_metabs.append(reac_or_metab.id[2:])

    reac_or_metab_list_df = pd.Series(list_of_reacs_or_metabs)
    reac_or_metab_list_df = pd.DataFrame(reac_or_metab_list_df, columns=[col_name])

    return reac_or_metab_list_df


def create_gpr_from_locus_tag(model: libModel, locus_tag: str, email: str) -> tuple[GeneProduct, libModel]:
    """Creates GeneProduct in the given model

    Args:
        - model (libModel): Model loaded with libSBML
        - locus_tag (str): NCBI compatible locus_tag
        - email (str): User Email to access the NCBI Entrez database

    Returns:
        tuple: libSBML GeneProduct (1) & libSBML model (2)
            (1) GeneProduct: Created gene product
            (2) libModel: Model containing the created gene product
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


def create_gp(model: libModel, model_id: str, name: str, locus_tag: str, protein_id: str) -> tuple[GeneProduct, libModel]:
    """Creates GeneProduct in the given model

    Args:
        - model (libModel): Model loaded with libSBML
        - model_id (str): ID identical to ID that CarveMe adds from the NCBI FASTA input file
        - name (str): Name of the GeneProduct
        - locus_tag (str): Genome-specific locus tag used as label in the model
        - protein_id (str): NCBI Protein/RefSeq ID

    Returns:
        tuple: libSBML GeneProduct (1) & libSBML model (2)
            (1) GeneProduct: Created gene product
            (2) libModel: Model containing the created gene product
    """
    id_db = None
    gp = model.getPlugin(0).createGeneProduct()
    gp.setId(model_id)
    gp.setName(name)
    gp.setLabel(locus_tag)
    gp.setSBOTerm('SBO:0000243')
    gp.setMetaId(f'meta_{model_id}')
    if re.fullmatch('^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', protein_id, re.IGNORECASE):
        id_db = 'REFSEQ'
    elif re.fullmatch('^(\w+\d+(\.\d+)?)|(NP_\d+)$', protein_id, re.IGNORECASE): id_db = 'NCBI'
    if id_db: add_cv_term_genes(protein_id, id_db, gp)
    return gp, model


def create_species(
    model: libModel, metabolite_id: str, name: str, compartment_id: str, charge: int, chem_formula: str
                   ) -> tuple[Species, libModel]:
    """Creates Species/Metabolite in the given model

    Args:
        - model (libModel): Model loaded with libSBML
        - metabolite_id (str): Metabolite ID within model (If model from CarveMe, preferable a BiGG ID)
        - name (str): Name of the metabolite
        - compartment_id (str): ID of the compartment where metabolite resides
        - charge (int): Charge for the metabolite
        - chem_formula (str): Chemical formula for the metabolite

    Returns:
        tuple: libSBML Species (1) & libSBML model (2)
        (1) Species: Created species/metabolite
        (2) libModel: Model containing the created metabolite
    """
    metabolite = model.createSpecies()
    metabolite.setId(f'M_{metabolite_id}')
    if name: metabolite.setName(name)
    metabolite.setMetaId(f'meta_M_{metabolite_id}')
    metabolite.setSBOTerm('SBO:0000247')
    metabolite.setInitialAmount(float('NaN'))
    metabolite.setHasOnlySubstanceUnits(True)
    metabolite.setBoundaryCondition(False)
    metabolite.setConstant(False)
    metabolite.setCompartment(compartment_id)
    metabolite.getPlugin(0).setCharge(charge)
    metabolite.getPlugin(0).setChemicalFormula(chem_formula)
    add_cv_term_metabolites(metabolite_id[:-2], 'BIGG', metabolite)
    return metabolite, model


def get_reversible(fluxes: dict[str: str]) -> bool:
    """Infer if reaction is reversible from flux bounds
    
    Args:
        - fluxes (dict): Dictionary containing the keys 'lower_bound' & 'upper_bound' 
                        with values in ['cobra_default_lb', 'cobra_0_bound', 'cobra_default_ub']
    Returns:
        bool: True if reversible else False
    """
    return (fluxes['lower_bound'] == 'cobra_default_lb') and (fluxes['upper_bound'] == 'cobra_default_ub')


def create_reaction(
    model: libModel, reaction_id: str, name:str, reactants: dict[str: int], products: dict[str: int], 
    fluxes: dict[str: str], reversible: bool=None, fast: bool=None, compartment: str=None, sbo: str=None, 
    genes: Union[str, list[str]]=None
    ) -> tuple[Reaction, libModel]:
    """Creates new reaction in the given model

    Args:
        - model (libModel): Model loaded with libSBML
        - reaction_id (str): BiGG ID of the reaction to create
        - name (str): Human readable name of the reaction
        - reactants (dict): Metabolites as keys and their stoichiometry as values
        - products (dict): Metabolites as keys and their stoichiometry as values
        - fluxes (dict): Dictionary with lower_bound and upper_bound as keys
        - reversible (bool): True/False for the reaction
        - fast (bool): True/False for the reaction
        - compartment (str): BiGG compartment ID of the reaction (if available)
        - sbo (str): SBO term of the reaction
        - genes (str|list): List of genes belonging to reaction

    Returns:
        tuple: libSBML reaction (1) & libSBML model (2)
            (1) Reaction: Created reaction 
            (2) libModel: Model containing the created reaction
    """
    reaction = model.createReaction()
    reaction.setId('R_' + reaction_id)
    if name: reaction.setName(name)
    reaction.setMetaId('meta_R_' + reaction_id)
    sbo = sbo if sbo else 'SBO:0000167'  # SBO term for biochemical or transport reaction
    reaction.setSBOTerm(sbo)
    fast = fast if fast else False
    reaction.setFast(fast)
    if compartment: reaction.setCompartment(compartment)  # Set compartment for reaction if available
    reversible = reversible if reversible else get_reversible(fluxes)
    reaction.setReversible(reversible)
    if genes:
            if genes == 'G_spontaneous':
                reaction.getPlugin(0).createGeneProductAssociation().createGeneProductRef().setGeneProduct(gene)
            elif len(genes) == 1:
                reaction.getPlugin(0).createGeneProductAssociation().createGeneProductRef().setGeneProduct(genes[0])
            else:
                gp_ass_or = reaction.getPlugin(0).createGeneProductAssociation().createOr()
                for gene in genes:
                    # Set GeneProductReferences if available
                    gp_ass_or.createGeneProductRef().setGeneProduct(gene)
    for metab, stoich in reactants.items(): #reactants as dict with metab:stoich
        reaction.addReactant(model.getSpecies('M_' + metab), stoich)
    for metab, stoich in products.items(): #reactants as dict with metab:stoich
        reaction.addProduct(model.getSpecies('M_' + metab), stoich)
    reaction.getPlugin(0).setLowerFluxBound(fluxes['lower_bound'])
    reaction.getPlugin(0).setUpperFluxBound(fluxes['upper_bound'])
    add_cv_term_reactions(reaction_id, 'BIGG', reaction)
    return reaction, model
 
