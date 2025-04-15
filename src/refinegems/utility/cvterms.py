#!/usr/bin/env python
"""Helper module to work with annotations (CVTerms)

Stores dictionaries which hold information the identifiers.org syntax, has functions to add CVTerms to different entities and parse CVTerms.
"""

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################

import cobra
import logging

from bioregistry import get_identifiers_org_iri, parse_iri
from libsbml import (
    BIOLOGICAL_QUALIFIER,
    BQB_IS,
    BQB_OCCURS_IN,
    BQB_IS_HOMOLOG_TO,
    MODEL_QUALIFIER,
    BQM_IS_DESCRIBED_BY,
)
from libsbml import Unit, CVTerm, Species, Reaction, GeneProduct, Group, SBase
from typing import Union

from ..developement.decorators import debug

################################################################################
# variables
################################################################################

DB2PREFIX_METABS = {
    "BIGG": "bigg.metabolite",
    "BIOCYC": "biocyc",
    "BRENDA": "brenda",
    "CHEBI": "chebi",
    "HMDB": "hmdb",
    "HUMAN METABOLOME DATABASE": "hmdb",
    "INCHI": "inchi",
    "INCHI-KEY": "inchikey",
    "KEGG": "kegg.compound",
    "KEGG COMPOUND": "kegg.compound",
    "METACYC": "metacyc.compound",
    "MNXREF": "metanetx.chemical",
    "METANETX": "metanetx.chemical",
    "PUBCHEM": "pubchem.compound",
    "REACTOME": "reactome",
    "SEED": "seed.compound",
    #"UPA": "unipathway.compound", #closed due to financial issues
    #"UNIPATHWAY COMPOUND": "unipathway.compound"
    "VMH": "vmhmetabolite",
}  #: :meta hide-value:

PREFIX2DB_METABS = {v: k for (k, v) in DB2PREFIX_METABS.items()}  #: :meta hide-value:

DB2PREFIX_REACS = {
    "BIGG": "bigg.reaction",
    "BIOCYC": "biocyc",
    "BRENDA": "brenda",
    "EC": "ec",
    "HMDB": "hmdb",
    "HUMAN METABOLOME DATABASE": "hmdb",
    "KEGG": "kegg.reaction",
    "METACYC": "metacyc.reaction",
    "MNXREF": "metanetx.reaction",
    "METANETX": "metanetx.reaction",
    "REACTOME": "reactome",
    "RHEA": "rhea",
    "SEED": "seed.reaction",
    #"UPA": "unipathway.reaction:", #closed due to financial issues
    #"UNIPATHWAY REACTION": "unipathway.reaction:"
    "VMH": "vmhreaction",
}  #: :meta hide-value:

PREFIX2DB_REACS = {v: k for (k, v) in DB2PREFIX_REACS.items()}  #: :meta hide-value:

DB2PREFIX_GENES = {
    "KEGG": "kegg.genes",
    "NCBI": "ncbiprotein",
    "NCBIGENE": "ncbigene",
    "REFSEQ": "refseq",
    "UNIPROT": "uniprot",
}  #: :meta hide-value:

PREFIX2DB_GENES = {v: k for (k, v) in DB2PREFIX_GENES.items()}  #: :meta hide-value:

DB2PREFIX_PATHWAYS = {"KEGG": "kegg.pathway"}  #: :meta hide-value:

PREFIX2DB_PATHWAYS = {
    v: k for (k, v) in DB2PREFIX_PATHWAYS.items()
}  #: :meta hide-value:

MIRIAM = "https://identifiers.org/"  #: :meta hide-value:
OLD_MIRIAM = "http://identifiers.org/"  #: :meta hide-value:

################################################################################
# functions
################################################################################

# cobra
# -----


def _add_annotations_from_dict_cobra(
    references: dict, entity: Union[cobra.Reaction, cobra.Metabolite, cobra.Model]
) -> None:
    """Given a dictionary and a cobra object, add the former as annotations to the latter.
    The keys of the dictionary are used as the annotation labels, the values as the values.
    If the keys are already in the entity, the values will be combined (union).

    Args:
        - references (dict):
            The dictionary with the references to add the entity.
        - entity (cobra.Reaction | cobra.Metabolite | cobra.Model):
            The entity to add annotations to.
    """
    # add additional references from the parameter
    for db, idlist in references.items():
        if not isinstance(idlist, list):
            idlist = [idlist]
        if db in entity.annotation.keys():
            if isinstance(entity.annotation[db], list):
                entity.annotation[db] = list(set(entity.annotation[db] + idlist))
            else:
                entity.annotation[db] = list(set(list(entity.annotation[db]) + idlist))
        else:
            entity.annotation[db] = idlist


# libsbml
# -------


def add_cv_term_units(unit_id: str, unit: Unit, relation: int):
    """Adds CVTerm to a unit

    Args:
        - unit_id (str):
            ID to add as URI to annotation
        - unit (Unit):
            Unit to add CVTerm to
        - relation (int):
            Provides model qualifier to be added
    """
    db_id = "pubmed" if relation == BQM_IS_DESCRIBED_BY else "uo"
    resource = get_identifiers_org_iri(db_id, unit_id)
    if resource:
        cv = CVTerm()
        cv.setQualifierType(MODEL_QUALIFIER)
        cv.setModelQualifierType(relation)
        cv.addResource(resource)
        unit.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {unit} with relation {relation} and {db_id}:{unit_id}."
        )


def add_cv_term_metabolites(entry: str, db_id: str, metab: Species):
    """Adds CVTerm to a metabolite

    Args:
        - entry (str):
            Id to add as annotation
        - db_id (str):
            Database to which entry belongs. Must be in DB2PREFIX_METABS.keys().
        - metab (Species):
            Metabolite to add CVTerm to
    """
    resource = get_identifiers_org_iri(DB2PREFIX_METABS.get(db_id), entry)
    if resource:
        if db_id == "HMDB" or db_id == "HUMAN METABOLOME DATABASE":
            if entry[:4] == "HMDB":
                entry = entry[4:]
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(resource)
        metab.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {metab} with {db_id} and {DB2PREFIX_METABS.get(db_id)}:{entry}."
        )


def add_cv_term_reactions(entry: str, db_id: str, reac: Reaction):
    """Adds CVTerm to a reaction

    Args:
        - entry (str):
            Id to add as annotation
        - db_id (str):
            Database to which entry belongs. Must be in DB2PREFIX_REACS.keys().
        - reac (Reaction):
            Reaction to add CVTerm to
    """
    resource = get_identifiers_org_iri(DB2PREFIX_REACS.get(db_id), entry)
    if resource:
        if db_id == "HMDB" or db_id == "HUMAN METABOLOME DATABASE":
            if entry[:4] == "HMDB":
                entry = entry[4:]
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(resource)
        reac.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {reac} with {db_id} and {DB2PREFIX_REACS.get(db_id)}:{entry}."
        )


def add_cv_term_genes(
    entry: str, db_id: str, gene: GeneProduct, lab_strain: bool = False
):
    """Adds CVTerm to a gene

    Args:
        - entry (str):
            Id to add as annotation.
        - db_id (str):
            Database to which entry belongs. Must be in DB2PREFIX_GENES.keys().
        - gene (GeneProduct):
            Gene to add CVTerm to.
        - lab_strain (bool, optional):
            For locally sequenced strains the qualifiers are always HOMOLOG_TO. Defaults to False.
    """
    resource = get_identifiers_org_iri(DB2PREFIX_GENES.get(db_id), entry)
    if resource:
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        if lab_strain:
            cv.setBiologicalQualifierType(BQB_IS_HOMOLOG_TO)
        else:
            cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(resource)
        gene.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {gene} with {db_id} and {DB2PREFIX_GENES.get(db_id)}:{entry}."
        )


def add_cv_term_pathways(entry: str, db_id: str, path: Group):
    """Add CVTerm to a groups pathway

    Args:
        - entry (str):
            Id to add as annotation
        - db_id (str):
            Database to which entry belongs. Must be in DB2PREFIX_PATHWAYS.keys().
        - path (Group):
            Pathway to add CVTerm to
    """
    resource = get_identifiers_org_iri(DB2PREFIX_PATHWAYS.get(db_id), entry)
    if resource:
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(resource)
        path.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {path} with {db_id} and {DB2PREFIX_PATHWAYS.get(db_id)}:{entry}."
        )


def add_cv_term_pathways_to_entity(entry: str, db_id: str, reac: Reaction):
    """Add CVTerm to a reaction as OCCURS IN pathway

    Args:
        - entry (str):
            Id to add as annotation
        - db_id (str):
            Database to which entry belongss
        - reac (Reaction):
            Reaction to add CVTerm to
    """
    resource = get_identifiers_org_iri(DB2PREFIX_PATHWAYS.get(db_id), entry)
    if resource:
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_OCCURS_IN)
        cv.addResource(resource)
        reac.addCVTerm(cv)
    else:
        logging.warning(
            f"No valid IRI could be formed for {reac} with {db_id} and {DB2PREFIX_PATHWAYS.get(db_id)}:{entry}."
        )

def get_id_from_cv_term(entity: SBase, db_id: str) -> list[str]:
    """Extract Id for a specific database from CVTerm

    Args:
        - entity (SBase):
            Species, Reaction, Gene, Pathway
        - db_id (str):
            Database of interest

    Returns:
        list[str]:
            Ids of entity belonging to db_id
    """
    def _extract_id_from_uri(uri: str) -> str:
        """Extracts the ID from the URI.

        Args:
            uri (str): 
                The URI to extract the ID from.

        Returns:
            str: 
                The extracted ID.
        """
        extracted_id = parse_iri(uri)[1]

        if not extracted_id:
            logging.info(f"Could not extract ID with bioregistry from URI: {uri}")
            extracted_id = uri.split('/')[-1]  # Fallback to splitting by "/"
            if ':' in extracted_id:
                extracted_id = extracted_id.split(':')[-1]

            if not extracted_id:
                logging.warning(f"Could not extract ID from URI: {uri} at all!")
                return None

        return extracted_id
        
    num_cvs = entity.getNumCVTerms()
    all_ids = []
    for i in range(0, num_cvs):
        ann_string = entity.getCVTerm(i)
        num_res = ann_string.getNumResources()
        ids = [
            _extract_id_from_uri(ann_string.getResourceURI(r))
            for r in range(0, num_res)
            if str(db_id) in ann_string.getResourceURI(r)
        ]
        all_ids.extend(ids)

    # Clean-up: Remove all Nones
    all_ids = [_ for _ in all_ids if _ is not None]
    if len(all_ids) == 0:
        logging.info(f'No URIs extracted for {db_id} database from {entity.getId()}')
    
    return all_ids


def generate_cvterm(qt, b_m_qt) -> CVTerm:
    """Generates a CVTerm with the provided qualifier & biological or model qualifier types

    Args:
        - qt (libSBML qualifier type):
            BIOLOGICAL_QUALIFIER or MODEL_QUALIFIER
        - b_m_qt (libSBML qualifier):
            BQM_IS, BQM_IS_HOMOLOG_TO, etc.

    Returns:
        CVTerm:
            With provided qualifier & biological or model qualifier types
    """
    cvterm = CVTerm()
    cvterm.setQualifierType(qt)

    if qt == BIOLOGICAL_QUALIFIER:
        cvterm.setBiologicalQualifierType(b_m_qt)
    else:
        cvterm.setModelQualifierType(b_m_qt)

    return cvterm


@debug
def print_cvterm(cvterm: CVTerm):
    """Debug function: Prints the URIs contained in the provided CVTerm along with the provided qualifier & biological/model qualifier types

    Args:
        cvterm (CVTerm):
            A libSBML CVTerm
    """
    if cvterm == None:
        logging.info("CVTerm currently empty!")
    else:
        current_b_m_qt = 0

        current_qt = cvterm.getQualifierType()

        if current_qt == BIOLOGICAL_QUALIFIER:
            current_b_m_qt = cvterm.getBiologicalQualifierType()
        elif current_qt == MODEL_QUALIFIER:
            current_b_m_qt = cvterm.getModelQualifierType()

        if cvterm.getNumResources() == 0:
            logging.info("No URIs present.")
        else:
            logging.info(f"Current CVTerm contains:  {cvterm.getResourceURI(0)}")

            for i in range(1, cvterm.getNumResources()):
                logging.info(f"                          {cvterm.getResourceURI(i)}")

        logging.info(
            f"Current CVTerm has QualifierType {current_qt} and Biological/ModelQualifierType {current_b_m_qt}."
        )
