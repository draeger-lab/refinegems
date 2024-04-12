#!/usr/bin/env python
""" Helper module to work with annotations (CVTerms)

Stores dictionaries which hold information the identifiers.org syntax, has functions to add CVTerms to different entities and parse CVTerms.
"""
import logging
from libsbml import BIOLOGICAL_QUALIFIER, BQB_IS, BQB_OCCURS_IN, BQB_IS_HOMOLOG_TO, MODEL_QUALIFIER, BQM_IS_DESCRIBED_BY, Unit, CVTerm, Species, Reaction, GeneProduct, Group, SBase

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel"

metabol_db_dict = {
                   'BIGG': 'bigg.metabolite:',
                   'BIOCYC': 'biocyc:META:',
                   'BioCyc': 'biocyc:META:',
                   'BRENDA': 'brenda:',
                   'CHEBI': 'CHEBI:',
                   'ChEBI': 'CHEBI:',
                   'HMDB': 'hmdb:HMDB',
                   'Human Metabolome Database': 'hmdb:HMDB',
                   'INCHI': 'inchi:',
                   'InChI': 'inchi:',
                   'InChI-Key': 'inchikey:',
                   'KEGG': 'kegg.compound:',
                   #'KEGG Compound': 'kegg.compound:',
                   'METACYC': 'metacyc.compound:',
                   'MXNREF': 'metanetx.chemical:',
                   'MetaNetX': 'metanetx.chemical:',
                   'PUBCHEM': 'pubchem.compound:',
                   'REACTOME': 'reactome:',
                   'Reactome': 'reactome:',
                   'SEED': 'seed.compound:',
                   #'UPA': 'unipathway.compound:', #closed due to financial issues
                   #'UniPathway Compound': 'unipathway.compound:'
                   'VMH': 'vmhmetabolite:'
                  }

reaction_db_dict = {
                    'BIGG': 'bigg.reaction:',
                    'BioCyc': 'biocyc:META:',
                    'BRENDA': 'brenda:',
                    'EC': 'ec-code:',
                    'HMDB': 'hmdb:HMDB',
                    'KEGG': 'kegg.reaction:',
                    'METACYC': 'metacyc.reaction:',
                    'MXNREF': 'metanetx.reaction:',
                    'MetaNetX': 'metanetx.reaction:',
                    'REACTOME': 'reactome:',
                    'Reactome': 'reactome:',
                    'RHEA': 'rhea:',
                    'SEED': 'seed.reaction:',
                    #'UPA': 'unipathway.reaction:',
                    #'UniPathway Reaction': 'unipathway.reaction:'
                    'VMH': 'vmhreaction:'
                   }

gene_db_dict = {
                'KEGG': 'kegg.genes:',
                'NCBI': 'ncbiprotein:',
                'REFSEQ': 'refseq:',
                'UNIPROT': 'uniprot:'
               }

pathway_db_dict = {'KEGG': 'kegg.pathway:'}

MIRIAM = 'https://identifiers.org/'
OLD_MIRIAM = 'http://identifiers.org/'


def add_cv_term_units(unit_id: str, unit: Unit, relation: int):
    """Adds CVTerm to a unit

    Args:
        - unit_id (str): ID to add as URI to annotation
        - unit (Unit): Unit to add CVTerm to
        - relation (int): Provides model qualifier to be added
    """
    cv = CVTerm()
    cv.setQualifierType(MODEL_QUALIFIER)
    cv.setModelQualifierType(relation)
    
    if relation == BQM_IS_DESCRIBED_BY:
        cv.addResource(f'https://identifiers.org/{unit_id}')
    else:
        cv.addResource(f'https://identifiers.org/UO:{unit_id}')
    unit.addCVTerm(cv)


def add_cv_term_metabolites(entry: str, db_id: str, metab: Species):
    """Adds CVTerm to a metabolite

    Args:
        - entry (str): Id to add as annotation
        - db_id (str): Database to which entry belongs. Must be in metabol_db_dict.keys().
        - metab (Species): Metabolite to add CVTerm to
    """
    if db_id == 'HMDB' or db_id == 'Human Metabolome Database':
        if entry[:4] == 'HMDB':
            entry = entry[4:]
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + metabol_db_dict[db_id] + entry)
    metab.addCVTerm(cv)
    metab.addCVTerm(cv)


def add_cv_term_reactions(entry: str, db_id: str, reac: Reaction):
    """Adds CVTerm to a reaction

    Args:
        - entry (str): Id to add as annotation
        - db_id (str): Database to which entry belongs. Must be in reaction_db_dict.keys().
        - reac (Reaction): Reaction to add CVTerm to
    """
    if db_id == 'HMDB' or db_id == 'Human Metabolome Database':
        if entry[:4] == 'HMDB':
            entry = entry[4:]
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource(
        'https://identifiers.org/' +
        reaction_db_dict[db_id] +
        entry)
    reac.addCVTerm(cv)


def add_cv_term_genes(entry: str, db_id: str, gene: GeneProduct, lab_strain: bool=False):
    """Adds CVTerm to a gene

    Args:
        - entry (str): Id to add as annotation
        - db_id (str): Database to which entry belongs. Must be in gene_db_dict.keys().
        - gene (GeneProduct): Gene to add CVTerm to
        - lab_strain (bool, optional): For locally sequenced strains the qualifiers are always HOMOLOG_TO. Defaults to False.
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    if lab_strain:
        cv.setBiologicalQualifierType(BQB_IS_HOMOLOG_TO)
    else:
        cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + gene_db_dict[db_id] + entry)
    gene.addCVTerm(cv)


def add_cv_term_pathways(entry: str, db_id: str, path: Group):
    """Add CVTerm to a groups pathway

    Args:
        - entry (str): Id to add as annotation
        - db_id (str): Database to which entry belongs. Must be in pathway_db_dict.keys().
        - path (Group): Pathway to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + pathway_db_dict[db_id] + entry)
    path.addCVTerm(cv)
    

def add_cv_term_pathways_to_entity(entry: str, db_id: str, reac: Reaction):
    """Add CVTerm to a reaction as OCCURS IN pathway

    Args:
        - entry (str): Id to add as annotation
        - db_id (str): Database to which entry belongss
        - reac (Reaction): Reaction to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_OCCURS_IN)
    cv.addResource('https://identifiers.org/' + pathway_db_dict[db_id] + entry)
    reac.addCVTerm(cv)


def get_id_from_cv_term(entity: SBase, db_id: str) -> list[str]:
    """Extract Id for a specific database from CVTerm

    Args:
        - entity (SBase): Species, Reaction, Gene, Pathway
        - db_id (str): Database of interest

    Returns:
        list[str]: Ids of entity belonging to db_id
    """
    num_cvs = entity.getNumCVTerms()
    all_ids = []
    for i in range(0, num_cvs):
        ann_string = entity.getCVTerm(i)
        num_res = ann_string.getNumResources()
        ids = [ann_string.getResourceURI(r).split(
            '/')[-1] for r in range(0, num_res) if str(db_id) in ann_string.getResourceURI(r)]
        ids = [id_string.split(':')[-1] for id_string in ids if ':' in  id_string]
        all_ids.extend(ids)

    return all_ids


def generate_cvterm(qt, b_m_qt) -> CVTerm:
    """Generates a CVTerm with the provided qualifier & biological or model qualifier types

    Args:
        - qt (libSBML qualifier type): BIOLOGICAL_QUALIFIER or MODEL_QUALIFIER
        - b_m_qt (libSBML qualifier): BQM_IS, BQM_IS_HOMOLOG_TO, etc.

    Returns:
        CVTerm: With provided qualifier & biological or model qualifier types 
    """
    cvterm = CVTerm()
    cvterm.setQualifierType(qt)
    
    if qt == BIOLOGICAL_QUALIFIER:
        cvterm.setBiologicalQualifierType(b_m_qt)
    else:
        cvterm.setModelQualifierType(b_m_qt)
    
    return cvterm 


def print_cvterm(cvterm: CVTerm):
    """Debug function: Prints the URIs contained in the provided CVTerm along with the provided qualifier & biological/model qualifier types

    Args:
        cvterm (CVTerm): A libSBML CVTerm
    """
    if cvterm == None:
        logging.info('CVTerm currently empty!')
    else:
        current_b_m_qt = 0
    
        current_qt = cvterm.getQualifierType()
                
        if current_qt == BIOLOGICAL_QUALIFIER:
            current_b_m_qt = cvterm.getBiologicalQualifierType()
        elif current_qt == MODEL_QUALIFIER:
            current_b_m_qt = cvterm.getModelQualifierType()
        
        if cvterm.getNumResources() == 0:
            logging.info('No URIs present.')
        else:
            logging.info(f'Current CVTerm contains:  {cvterm.getResourceURI(0)}')
    
            for i in range(1, cvterm.getNumResources()):
                logging.info(f'                          {cvterm.getResourceURI(i)}')
                
        logging.info(f'Current CVTerm has QualifierType {current_qt} and Biological/ModelQualifierType {current_b_m_qt}.')