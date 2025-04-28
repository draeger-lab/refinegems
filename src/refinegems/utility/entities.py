#!/usr/bin/env python
"""Collection of functions to access, handle and manipulate different entities of
COBRApy and libsbml models.
"""

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################
import cobra
from cobra.io.sbml import _f_gene_rev
import json
import logging
import pandas as pd
import re
import requests
import urllib
import warnings

from Bio import Entrez
from Bio.KEGG import REST, Compound
from cobra.io.sbml import _f_gene
from datetime import date
from functools import reduce
from libsbml import Model as libModel
from libsbml import (
    Species,
    Reaction,
    FbcOr,
    FbcAnd,
    GeneProductRef,
    Unit,
    UnitDefinition,
    ListOfUnitDefinitions,
)
from libsbml import UNIT_KIND_MOLE, UNIT_KIND_GRAM, UNIT_KIND_LITRE, UNIT_KIND_SECOND
from libsbml import BQM_IS, BQM_IS_DERIVED_FROM, BQM_IS_DESCRIBED_BY
from pathlib import Path

from random import choice
from string import ascii_uppercase, digits
from typing import Union, Literal
from tqdm import tqdm

tqdm.pandas()
pd.options.mode.chained_assignment = None  # suppresses the pandas SettingWithCopyWarning; comment out before developing!!

from .cvterms import (
    add_cv_term_genes,
    add_cv_term_metabolites,
    add_cv_term_reactions,
    add_cv_term_units,
    _add_annotations_from_dict_cobra,
)
from .db_access import (
    kegg_reaction_parser,
    _add_annotations_from_bigg_reac_row,
    get_BiGG_metabs_annot_via_dbid,
    add_annotations_from_BiGG_metabs,
    _search_ncbi_for_gp,
)
from .io import load_a_table_from_database, parse_gff_for_cds
from .util import (
    COMP_MAPPING,
    DB2REGEX,
    VALID_COMPARTMENTS,
    MIN_GROWTH_THRESHOLD,
    is_stoichiometric_factor,
)
from ..developement.decorators import *

################################################################################
# variables
################################################################################

REF_COL_GF_GENE_MAP = {"UniProt": "UNIPROT", "GeneID": "NCBIGENE"}  #: :meta:

################################################################################
# functions
################################################################################

# ++++++++++++++++++++++++++++++++++++++++
# COBRApy - models
# ++++++++++++++++++++++++++++++++++++++++

# handling compartments
# ---------------------


def are_compartment_names_valid(model: cobra.Model) -> bool:
    """Check if compartment names of model are considered valid based on
    :py:const:`~refinegems.utility.util.VALID_COMPARTMENTS`.

    Args:
        - model (cobra.Model):
            The model, loaded with COBRApy.

    Returns:
        bool:
            True, if valid, else false.
    """

    for c in model.compartments.keys():
        if c not in VALID_COMPARTMENTS.keys():
            return False

    return True


def resolve_compartment_names(model: cobra.Model):
    """Resolves compartment naming problems.

    Args:
        - model (cobra.Model):
            A COBRApy model object.

    Raises:
        - KeyError: Unknown compartment raises an error to add it to the mapping. Important for developers.
    """

    # check if compartment names are valid
    if not are_compartment_names_valid(model):

        # check if mapping is possible
        if set(model.compartments.keys()).issubset(set(COMP_MAPPING.keys())):
            # for each metabolite rename the compartment
            for metabolite in model.metabolites:
                metabolite.compartment = COMP_MAPPING[metabolite.compartment]
            # add whole descriptions of the compartments to the model
            # note:
            #    only compartments IN the model will be added
            model.compartments = VALID_COMPARTMENTS
            if "uc" in model.compartments.keys():
                logging.warning("Unknown compartment(s) detected and (re)named 'uc'. Simulation results might be affected.")

        else:
            raise KeyError(
                f"Unknown compartment {[_ for _ in model.compartments if _ not in COMP_MAPPING.keys()]} detected. Please contact developers or change compartment names manually."
            )


# handling cobra entities (features)
# ----------------------------------


def reaction_equation_to_dict(eq: str, model: cobra.Model) -> dict:
    """Parses a reaction equation string to dictionary

    Args:
        - eq (str):
            Equation of a reaction
        - model (cobra.Model):
            Model loaded with COBRApy

    Returns:
        dict:
            Metabolite Ids as keys and their coefficients as values (negative = reactants, positive = products)
    """
    # from Alina Renz
    eq = eq.split(" ")
    eq_matrix = {}
    are_products = False
    coeff = 1
    for i, part in enumerate(eq):
        if part == "-->":
            are_products = True
            continue
        if part == "+":
            continue
        if part == "2":
            coeff = 2
            continue
        if are_products:
            eq_matrix[model.metabolites.get_by_id(part)] = 1 * coeff
            coeff = 1
        else:
            eq_matrix[model.metabolites.get_by_id(part)] = -1 * coeff
            coeff = 1
    return eq_matrix


def get_reaction_annotation_dict(
    model: cobra.Model, db: Literal["KEGG", "BiGG"]
) -> dict:
    """Create a dictionary of a model's reaction IDs and a chosen database ID as
    saved in the annotations of the model.

    The database ID can be choosen based on the strings for the namespace options
    in other functions.

    Args:
        - model (cobra.Model):
            A model loaded with COBRApy.
        - db (Literal['KEGG','BiGG']):
            The string denoting the database to map to.

    Raises:
        - ValueError: Unknown database string for paramezer db

    Returns:
        dict:
            The mapping of the reaction IDs to the database IDs found in the annotations
    """

    react_dict = {}

    match db:
        case "KEGG":
            db_string = "kegg.reaction"
        case "BiGG":
            db_string = "bigg.reaction"
        case _:
            mes = f"Unknown database string for parameter db: {db}"
            raise ValueError(mes)

    for r in model.reactions:
        if db_string in r.annotation.keys():
            react_dict[r.id] = r.annotation[db_string]
        else:
            react_dict[r.id] = "-"

    return react_dict


def create_random_id(
    model: cobra.Model, entity_type: Literal["reac", "meta"] = "reac", prefix: str = ""
) -> str:
    """Generate a unique, random ID for a model entity for a model.

    Args:
        - model (cobra.Model):
            A model loaded with COBRApy.
        - entity_type (Literal['reac','meta'], optional):
            Type of model entity.
            Can be 'reac' for Reaction or 'meta' for Metabolite.
            Defaults to 'reac'.
        - prefix (str, optional):
            Prefix to set for the randomised part.
            Useful to identify the random IDs later on.
            Defaults to ''.

    Raises:
        - ValueError: Unknown entity_type

    Returns:
        str:
            The generate new and unique ID.
    """

    match entity_type:
        case "reac":
            all_ids = [_.id for _ in model.reactions]
        case "meta":
            all_ids = [_.id for _ in model.metabolites]
        case _:
            mes = f"Unkown entity_type: {entity_type}"
            raise ValueError(mes)

    prefix = f"{prefix}{entity_type}"
    var = "".join(choice(ascii_uppercase + digits) for i in range(6))
    label = prefix + var
    j = 6

    while True:

        for i in range(36**6):  # make sure it does not run endlessly
            if label in all_ids:
                label = prefix + "".join(
                    choice(ascii_uppercase + digits) for x in range(j)
                )
            else:
                return label

        j = j + 1


def match_id_to_namespace(
    model_entity: Union[cobra.Reaction, cobra.Metabolite], namespace: Literal["BiGG"]
) -> None:
    """Based on a given namespace, change the ID of a given model entity to the set namespace.

    Currently working namespaces:

        - BiGG

    Args:
        - model_entity (cobra.Reaction, cobra.Metabolite]):
            The model entity.
            Can be either a cobra.Reaction or cobra.Metabolite object.
        - namespace (Literal['BiGG']):
            The chosen namespace.

    Raises:
        - ValueError: Unknown input for namespace
        - TypeError: Unknown type for model_entity
    """

    match model_entity:

        # Reaction
        # --------
        case cobra.Reaction():
            match namespace:

                case "BiGG":
                    if "bigg.reaction" in model_entity.annotation.keys():
                        # currently takes first entry if annotation is list
                        model_entity.id = (
                            model_entity.annotation["bigg.reaction"]
                            if isinstance(model_entity.annotation["bigg.reaction"], str)
                            else model_entity.annotation["bigg.reaction"][0]
                        )

                case _:
                    mes = f"Unknown input for namespace: {namespace}"
                    raise ValueError(mes)

        # Metabolite
        # ----------
        case cobra.Metabolite():
            match namespace:

                case "BiGG":
                    if "bigg.metabolite" in model_entity.annotation.keys():
                        # currently takes first entry if annotation is list
                        model_entity.id = (
                            model_entity.annotation["bigg.metabolite"]
                            + "_"
                            + model_entity.compartment
                            if isinstance(
                                model_entity.annotation["bigg.metabolite"], str
                            )
                            else model_entity.annotation["bigg.metabolite"][0]
                        )

                case _:
                    mes = f"Unknown input for namespace: {namespace}"
                    raise ValueError(mes)
        # Error
        # -----
        case _:
            mes = f"Unknown type for model_entity: {type(model_entity)}"
            raise TypeError(mes)


def isreaction_complete(
    reac: cobra.Reaction,
    name_check: bool = False,
    formula_check: Literal["none", "existence", "wildcard", "strict"] = "existence",
    exclude_dna: bool = True,
    exclude_rna: bool = True,
) -> bool:
    """Check, if a reaction object can be considered a complete reaction.
    The parameters can set the strictness of the checking.
    Useful for checking, if a reaction can be added to a model
    or if it might break it. For example, a missing model ID in either reaction or
    metabolites returns false.



    Args:
        - reac (cobra.Reaction):
            The reaction object (COBRApy) to test.
        - name_check (bool, optional):
            Option to force reaction and metabolites to have the name attribute set.
            Defaults to False.
        - formula_check (Literal['none','existence','wildcard','strict'], optional):
            Option to check the formula. 'none' disables the check,
            'existence' tests, if a formula is set, 'wildcard' additionally checks for wild cards
            (returns false if one found) and 'strict' also checks for the rest symbol 'R'.
            Defaults to 'existence'.
        - exclude_dna (bool, optional):
            Option to set DNA reactions to invalid. Defaults to True.
        - exclude_rna (bool, optional):
            Option to set RNA reaction to invalid. Defaults to True.

    Returns:
        bool:
            True, if the checks are passed successfully, else False.
    """

    # check reaction features
    # -----------------------
    # check name
    if name_check:
        if reac.id == "" or pd.isnull(reac.name):
            return False
    # check for RNA/DNA
    if exclude_dna and "DNA" in reac.name:
        return False
    if exclude_rna and "RNA" in reac.name:
        return False

    # check metabolites features
    # --------------------------
    for m in reac.metabolites:
        # check id
        if m.id == "" or pd.isnull(m.id):
            return False
        # check name
        if name_check and (m.name == "" or pd.isnull(m.name)):
            return False
        # check formula
        match formula_check:
            # case 1: no formula checking
            case "none":
                pass
            # case 2: check, if formula is set
            case "existence":
                if m.formula == "" or pd.isnull(m.formula):
                    return False
            # case 3: check, if formula is set and only contains alphanumerical characters
            case "wildcard":
                if m.formula == "" or pd.isnull(m.formula) or not m.formula.isalnum():
                    return False
            # case 4: check, if formula is set, constains only alphanumerical characters and no rest R
            case "strict":
                if (
                    m.formula == ""
                    or pd.isnull(m.formula)
                    or not m.formula.isalnum()
                    or bool(re.search(r"R(?![a-z])", "CH2HOROH"))
                ):
                    return False
            # default case: Warning + no formula check
            case _:
                mes = f"Unknown options for formula_check: {formula_check}\nChecking the metabolite formula will be skipped."
                warnings.warn(mes, UserWarning)

    return True


# adding metabolites to cobra models
# ----------------------------------


@template
def build_metabolite_xxx(
    id: str, model: cobra.Model, namespace: str, compartment: str, idprefix: str
) -> cobra.Metabolite:
    """Template function for building a cobra.Metabolite.

    .. note::

        This is a template function for developers. It cannot be executed.

    Args:
        - id (str):
            _description_
        - model (cobra.Model):
            _description_
        - namespace (str):
            _description_
        - compartment (str):
            _description_
        - idprefix (str):
            _description_

    Returns:
        cobra.Metabolite:
            _description_
    """
    # change xxx to database or way the metabolite will be reconstructed with
    # check if id in model
    # get information via id
    # collection formation in a new metabolite object
    # add more annotations from other databases
    # adjust namespace
    # check model again for new namespace
    pass


# originally from SPECIMEN
def build_metabolite_mnx(
    id: str,
    model: cobra.Model,
    namespace: str = "BiGG",
    compartment: str = "c",
    idprefix: str = "refineGEMs",
) -> Union[cobra.Metabolite, None]:
    """Build a cobra.Metabolite object from a MetaNetX ID.
    This function will NOT directly add the metabolite to the model,
    if the contruction is successful.

    Args:
        - id (str):
            A MetaNetX ID of a metabolite.
        - model (cobra.Model):
            The model, the metabolite will be build for.
        - namespace (str, optional):
            Name to use for the model ID. If namespace cannot be matched,
            will use a random ID.
            Defaults to 'BiGG'.
        - compartment (str, optional):
            Compartment of the metabolite. Defaults to 'c'.
        - idprefix (str, optional):
            Prefix for the random ID. Defaults to 'refineGEMs'.

    Returns:
        1. Case construction successful or match found in model:
                cobra.Metabolite:
                    The metabolite object.

        2. Case construction failed:
                None:
                    Nothing to return.
    """

    # fast check if compound already in model
    # ------------------------------------------
    # step 1: check if MetaNetX ID in model
    matches = [
        x.id
        for x in model.metabolites
        if "metanetx.chemical" in x.annotation
        and x.annotation["metanetx.chemical"] == id
        and x.compartment == compartment
    ]

    # step 2: if yes, retrieve metabolite from model
    # case 1: multiple matches found
    if len(matches) > 0:
        if len(matches) > 1:
            # currently, just the first match is taken
            match = model.metabolites.get_by_id(matches[0])
        #  case 2: only one match found
        else:
            match = model.metabolites.get_by_id(matches[0])

        # step 3: add metabolite
        return match

    # if not, create new metabolite
    # -----------------------------
    metabolite_prop = load_a_table_from_database(
        f"SELECT * FROM mnx_chem_prop WHERE id = '{id}'"
    )
    metabolite_anno = load_a_table_from_database(
        f"SELECT * FROM mnx_chem_xref WHERE id = '{id}'"
    )
    if len(metabolite_prop) == 0:  # cannot construct metabolite
        return None
    else:

        # step 1: create a random metabolite ID
        new_metabolite = cobra.Metabolite(create_random_id(model, "meta", idprefix))

        # step 2: add features
        # --------------------
        new_metabolite.formula = metabolite_prop["formula"].iloc[0]
        new_metabolite.name = metabolite_prop["name"].iloc[0]
        new_metabolite.charge = metabolite_prop["charge"].iloc[0]
        new_metabolite.compartment = compartment

        # step 3: add notes
        # -----------------
        new_metabolite.notes["created with"] = "refineGEMs based on metanetx.chemical"

        # step 4: add annotations
        # -----------------------
        # add SBOTerm
        new_metabolite.annotation["sbo"] = "SBO:0000247"

        # add information directly available from the mnx_chem_prop table
        new_metabolite.annotation["metanetx.chemical"] = [metabolite_prop["id"].iloc[0]]
        if not pd.isnull(metabolite_prop["InChIKey"].iloc[0]):
            new_metabolite.annotation["inchikey"] = (
                metabolite_prop["InChIKey"].iloc[0].split("=")[1]
            )

        # get more annotation from the mnx_chem_xref table
        for db in [
            "kegg.compound",
            "metacyc.compound",
            "seed.compound",
            "bigg.metabolite",
            "chebi",
        ]:
            db_matches = metabolite_anno[metabolite_anno["source"].str.contains(db)]
            if len(db_matches) > 0:
                new_metabolite.annotation[db] = [
                    m.split(":", 1)[1] for m in db_matches["source"].tolist()
                ]

        # Cleanup BiGG annotations
        if not "bigg.metabolite" in new_metabolite.annotation.keys():
            # if no BiGG was found in MetaNetX, try reverse search in BiGG
            get_BiGG_metabs_annot_via_dbid(
                new_metabolite, id, "MetaNetX (MNX) Chemical", compartment
            )

            # add additional information from BiGG (if ID found)
            add_annotations_from_BiGG_metabs(new_metabolite)

        # step 5: change ID according to namespace
        # ----------------------------------------
        match_id_to_namespace(new_metabolite, namespace)

        # step 6: re-check existence of ID in model
        # -----------------------------------------
        if new_metabolite.id in [_.id for _ in model.metabolites]:
            return model.metabolites.get_by_id(new_metabolite.id)

    return new_metabolite


# originally from SPECIMEN
def build_metabolite_kegg(
    kegg_id: str,
    model: cobra.Model,
    namespace: Literal["BiGG"] = "BiGG",
    compartment: str = "c",
    idprefix="refineGEMs",
) -> Union[cobra.Metabolite, None]:
    """Build a cobra.Metabolite object from a KEGG ID.
    This function will NOT directly add the metabolite to the model,
    if the contruction is successful.

    Args:
        - kegg_id (str):
            A KEGG ID of a metabolite.
        - model (cobra.Model):
            The model, the metabolite will be build for.
        - namespace (str, optional):
            Name to use for the model ID. If namespace cannot be matched,
            will use a random ID.
            Defaults to 'BiGG'.
        - compartment (str, optional):
            Compartment of the metabolite. Defaults to 'c'.
        - idprefix (str, optional):
            Prefix for the random ID. Defaults to 'refineGEMs'.

    Returns:
        1. Case construction successful or match found in model:
                cobra.Metabolite:
                    The build metabolite object.

        2. Case construction failed:
                None:
                    Nothing to return.
    """

    # ---------------------------------------
    # fast check if compound already in model
    # ---------------------------------------
    # step 1: check via KEGG ID
    matches = [
        x.id
        for x in model.metabolites
        if (
            "kegg.compound" in x.annotation and x.annotation["kegg.compound"] == kegg_id
        )
    ]
    if len(matches) > 0:
        # step 2: model id --> metabolite object
        #  case 1: multiple matches found
        if len(matches) > 1:
            # Currently takes only the first match
            match = model.metabolites.get_by_id(matches[0])
        #  case 2: only one match found
        else:
            match = model.metabolites.get_by_id(matches[0])

        # step 3: add metabolite
        return match

    # -----------------------------
    # if not, create new metabolite
    # -----------------------------

    # step 1: retrieve KEGG entry for compound
    # ----------------------------------------
    try:
        kegg_handle = REST.kegg_get(kegg_id)
        kegg_record = [r for r in Compound.parse(kegg_handle)][0]
    except urllib.error.HTTPError:
        warnings.warn(f"HTTPError: {kegg_id}")
        return None
    except ConnectionResetError:
        warnings.warn(f"ConnectionResetError: {kegg_id}")
        return None
    except urllib.error.URLError:
        warnings.warn(f"URLError: {kegg_id}")
        return None

    # step 2: create a random metabolite ID
    # -------------------------------------
    new_metabolite = cobra.Metabolite(create_random_id(model, "meta", idprefix))

    # step 3: add features
    # --------------------
    # set name from KEGG and additionally use it as ID if there is none yet
    if isinstance(kegg_record.name, list) and len(kegg_record.name) > 0:
        new_metabolite.name = kegg_record.name[0]
    elif isinstance(kegg_record.name, str) and len(kegg_record.name) > 0:
        new_metabolite.name = kegg_record.name
    else:
        new_metabolite.name = ""
    # set compartment
    new_metabolite.compartment = compartment
    # set formula
    new_metabolite.formula = kegg_record.formula

    # step 4: add notes
    # -----------------
    new_metabolite.notes["created with"] = "refineGEMs based on KEGG.compound"

    # step 5: add annotations
    # -----------------------
    # add annotation from the KEGG entry
    new_metabolite.annotation["kegg.compound"] = kegg_id
    db_idtf = {"CAS": "cas", "PubChem": "pubchem.compound", "ChEBI": "chebi"}
    for db, ids in kegg_record.dblinks:
        if db in db_idtf:
            new_metabolite.annotation[db_idtf[db]] = ids

    # add SBOTerm
    new_metabolite.annotation["sbo"] = "SBO:0000247"

    # search for infos in MetaNetX
    mnx_info = load_a_table_from_database(
        f"SELECT * FROM mnx_chem_xref WHERE source = 'kegg.compound:{kegg_id}'",
        query=True,
    )
    # if matches have been found
    if len(mnx_info) > 0:
        mnx_ids = list(set(mnx_info["id"]))
        
        # in case of an ambiguous mapping, take first match
        if len(mnx_ids) > 1:
            logging.warning('MNX mapping ambiguously due to multiple matches: {mnx_ids}\n-> Using first entry for building metabolite.')
            mnx_ids = [mnx_ids[0]]
        
        # mapping is unambiguously (or made that way)
        mnx_info = load_a_table_from_database(
            f"SELECT * FROM mnx_chem_prop WHERE id = '{mnx_ids[0]}'", query=True
        )
        # add charge
        new_metabolite.charge = mnx_info["charge"].iloc[0]
        # add more annotations
        new_metabolite.annotation["metanetx.chemical"] = [mnx_info["id"].iloc[0]]
        if not pd.isnull(mnx_info["InChIKey"].iloc[0]):
            new_metabolite.annotation["inchikey"] = (
                mnx_info["InChIKey"].iloc[0].split("=")[1]
            )

        # get more annotation from the mnx_chem_xref table
        metabolite_anno = load_a_table_from_database(
            f'SELECT * FROM mnx_chem_xref WHERE id = \'{mnx_info["id"]}\''
        )
        for db in [
            "kegg.compound",
            "metacyc.compound",
            "seed.compound",
            "bigg.metabolite",
            "chebi",
        ]:
            db_matches = metabolite_anno[metabolite_anno["source"].str.contains(db)]
            if len(db_matches) > 0:
                mnx_tmp = [
                    m.split(":", 1)[1] for m in db_matches["source"].tolist()
                ]
                if db in new_metabolite.annotation.keys():
                    new_metabolite.annotation[db] = list(
                        set(mnx_tmp + new_metabolite.annotation[db])
                    )
                else:
                    new_metabolite.annotation[db] = mnx_tmp

    # Cleanup BiGG annotations
    if not "bigg.metabolite" in new_metabolite.annotation.keys():
        # if no BiGG was found in KEGG, try reverse search in BiGG
        get_BiGG_metabs_annot_via_dbid(new_metabolite, id, "KEGG Compound", compartment)

        # add additional information from BiGG (if ID found)
        add_annotations_from_BiGG_metabs(new_metabolite)

    # step 6: change ID according to namespace
    # ----------------------------------------
    match_id_to_namespace(new_metabolite, namespace)

    # step 7: re-check existence of ID in model
    # -----------------------------------------
    if new_metabolite.id in [_.id for _ in model.metabolites]:
        return model.metabolites.get_by_id(new_metabolite.id)

    return new_metabolite


def build_metabolite_bigg(
    id: str,
    model: cobra.Model,
    namespace: Literal["BiGG"] = "BiGG",
    idprefix: str = "refineGEMs",
) -> Union[cobra.Metabolite, None]:
    """Build a cobra.Metabolite object from a BiGG ID.
    This function will NOT directly add the metabolite to the model,
    if the contruction is successful.

    Args:
        - id (str):
            A BiGG ID of a metabolite.
        - model (cobra.Model):
            The model, the metabolite will be build for.
        - namespace (str, optional):
            Name to use for the model ID. If namespace cannot be matched,
            will use a random ID.
            Defaults to 'BiGG'.
        - compartment (str, optional):
            Compartment of the metabolite. Defaults to 'c'.
        - idprefix (str, optional):
            Prefix for the random ID. Defaults to 'refineGEMs'.

    Returns:
        1. Case construction successful or match found in model:
                cobra.Metabolite:
                    The build metabolite object.

        2. Case construction failed:
                None:
                    Nothing to return.
    """

    compartment = id.rsplit("_", 1)[1]
    # ------------------------------------------
    # fast check if compound already in model
    # ------------------------------------------
    # step 1: check if MetaNetX ID in model
    matches = [
        x.id
        for x in model.metabolites
        if "bigg.metabolite" in x.annotation
        and (
            x.annotation["bigg.metabolite"] == id
            or x.annotation["bigg.metabolite"] == id.rsplit("_", 1)[0]
        )
        and x.compartment == compartment
    ]
    # step 2: if yes, retrieve metabolite from model
    # case 1: multiple matches found
    if len(matches) > 0:
        if len(matches) > 1:
            # currently, just the first one is taken
            match = model.metabolites.get_by_id(matches[0])
        #  case 2: only one match found
        else:
            match = model.metabolites.get_by_id(matches[0])

        # step 3: add metabolite
        return match

    # -----------------------------
    # if not, create new metabolite
    # -----------------------------
    # get information from the database
    bigg_res = load_a_table_from_database(
        f"SELECT * FROM bigg_metabolites WHERE id = '{id}'", query=True
    )
    if len(bigg_res) > 0:
        bigg_res = bigg_res.iloc[0, :]
    else:
        return None  # not data = no recontruction
    # get information from MNX if ID available
    mnx_res = None
    if bigg_res["MetaNetX (MNX) Chemical"]:
        mnx_res = load_a_table_from_database(
            f'SELECT * FROM mnx_chem_prop WHERE id=\'{bigg_res["MetaNetX (MNX) Chemical"]}\'',
            query=True,
        )
        if len(mnx_res) > 0:
            mnx_res = mnx_res.iloc[0, :]
        else:
            mnx_res = None

    # step 1: create a random metabolite ID
    # -------------------------------------
    new_metabolite = cobra.Metabolite(create_random_id(model, "meta", idprefix))

    # step 2: add features
    # --------------------
    new_metabolite.name = bigg_res["name"]
    new_metabolite.compartment = compartment
    if mnx_res is not None:
        if mnx_res["charge"]:
            new_metabolite.charge = mnx_res["charge"]
        if mnx_res["formula"]:
            new_metabolite.formula = mnx_res["formula"]

    if not new_metabolite.formula or not new_metabolite.charge:
        try:
            bigg_fetch = json.loads(
                requests.get(
                    f'http://bigg.ucsd.edu/api/v2/universal/metabolites/{id.rsplit("_",1)[0]}'
                ).text
            )
            if "formulae" in bigg_fetch.keys() and not new_metabolite.formula:
                new_metabolite.formula = bigg_fetch["formulae"][0]
            if "charges" in bigg_fetch.keys() and new_metabolite.charge:
                new_metabolite.charge = bigg_fetch["charges"][0]
        except Exception as e:
            pass

    # step 3: add notes
    # -----------------
    new_metabolite.notes["created with"] = "refineGEMs based on BiGG"

    # step 4: add annotations
    # -----------------------
    # add SBOTerm
    new_metabolite.annotation["sbo"] = "SBO:0000247"
    # add infos from BiGG
    new_metabolite.annotation["bigg.metabolite"] = [id.removesuffix(f"_{compartment}")]
    add_annotations_from_BiGG_metabs(new_metabolite)
    # add annotations from MNX
    if mnx_res is not None:
        if mnx_res["InChI"] and "inchi" not in new_metabolite.annotation.keys():
            new_metabolite.annotation["inchi"] = mnx_res["InChI"]
        if mnx_res["InChIKey"] and "inchikey" not in new_metabolite.annotation.keys():
            new_metabolite.annotation["inchikey"] = mnx_res["InChIKey"]

    # step 5: change ID according to namespace
    # ----------------------------------------
    match_id_to_namespace(new_metabolite, namespace)

    # step 6: re-check existence of ID in model
    # -----------------------------------------
    if new_metabolite.id in [_.id for _ in model.metabolites]:
        return model.metabolites.get_by_id(new_metabolite.id)

    return new_metabolite


# adding reactions to cobra models
# --------------------------------


def parse_reac_str(
    equation: str, type: Literal["BiGG", "BioCyc", "MetaNetX", "KEGG"] = "MetaNetX"
) -> tuple[dict, dict, list, bool]:
    """Parse a reaction string.

    Args:
        - equation (str):
            The equation of a reaction as a string (as saved in the database).
        - type (Literal['BiGG','BioCyc','MetaNetX','KEGG'], optional):
            The name of the database the equation was taken from.
            Can be 'BiGG','BioCyc','MetaNetX','KEGG'.
            Defaults to 'MetaNetX'.

    Returns:
        tuple:
            Tuple of (1) dict, (2) dict, (3) list & (4) bool:

            1. Dictionary with the reactant IDs and their stoichiometric factors.
            2. Dictionary with the product IDs and their stoichiometric factors.
            3. List of compartment IDs or None, if they cannot be extract from the equation.
            4. True, if the reaction is reversible, else False.
    """

    products = {}
    reactants = {}
    compartments = list()
    is_product = False
    reversible = True

    match type:
        case "MetaNetX":
            for s in equation.split(" "):
                # switch from reactants to products
                if s == "=":
                    is_product = True
                # found stoichiometric factor
                elif s.isnumeric():
                    factor = float(s)
                # skip
                elif s == "+":
                    continue
                # found metabolite
                else:
                    # get information from MetaNetX
                    metabolite, compartment = s.split("@")
                    compartments.append(compartment)

                    if is_product:
                        products[metabolite] = factor
                    else:
                        reactants[metabolite] = factor

        case "BiGG":
            factor = 1.0  # BiGG does not use factor 1 in the equations
            for s in equation.split(" "):
                # found factor
                if is_stoichiometric_factor(s):
                    factor = float(s)
                # switch from reactants to products
                elif s == "-->":
                    is_product = True
                    reversible = False
                elif s == "<->":
                    is_product = True
                # skip
                elif s == "+":
                    continue
                # found metabolite
                else:
                    if s:
                        compartments.append(s.rsplit("_", 1)[1])
                        if is_product:
                            products[s] = factor
                        else:
                            reactants[s] = factor
                    factor = 1.0

        case "KEGG":
            compartments = None
            factor = 1.0
            for s in equation.split(" "):
                if s.isnumeric():
                    factor = float(s)
                elif s == "+":
                    continue
                elif s == "<=>":
                    is_product = True
                else:
                    if is_product:
                        products[s] = factor
                    else:
                        reactants[s] = factor
                    factor = 1.0

    return (reactants, products, compartments, reversible)


def validate_reaction_compartment_bigg(comps: list) -> Union[bool, Literal["exchange"]]:
    """Retrieves and validates the compatment(s) of a reaction from the BiGG namespace

    Args:
        - comps (list):
            List containing compartments from a BiGG reaction

    Returns:

        (1) Case ``compartment in VALID_COMPARTMENTS.keys()``
                bool|'exchange':
                    Either

                    - True if the provided compartments are valid
                    - 'exchange' if reaction in multiple compartments

        (2) Case not a valid compartment:
                bool:
                    'False' if one of the found compartments is not in :py:const:`~refinegems.utility.util.VALID_COMPARTMENTS`
    """

    contained_in_compartments = [
        (comp in VALID_COMPARTMENTS.keys()) for comp in comps
    ]  # Get True for correct compartment
    if not all(contained_in_compartments):  # At least one compartment not correct
        return False
    else:  # All compartments correct
        if (
            len(set(comps)) == 1
        ):  # Set of found compartments of reaction = 1: Reaction happens in one valid compartment
            return True
        else:  # Probably exchange reaction
            return "exchange"


@template
def build_reaction_xxx():
    """
    Extend the build function so, that all of them can take either the id or an equation
    as input for rebuilding the reaction (would also be beneficial for semi-manual curation)

    model:cobra.Model, id:str=None,
                       reac_str:str=None,
                       references:dict={},
                       idprefix:str='refineGEMs',
                       namespace:Literal['BiGG']='BiGG') -> Union[cobra.Reaction, None, list]:
    """
    pass


def build_reaction_mnx(
    model: cobra.Model,
    id: str,
    reac_str: str = None,
    references: dict = {},
    idprefix: str = "refineGEMs",
    namespace: Literal["BiGG"] = "BiGG",
) -> Union[cobra.Reaction, None, list]:
    """Construct a new reaction for a model from a MetaNetX reaction ID.
    This function will NOT add the reaction directly to the model, if the
    construction process is successful.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - id (str):
            A MetaNetX reaction ID.
        - reac_str (str, optional):
            The reaction equation string from the database.
            Defaults to None.
        - references (dict, optional):
            Additional annotations to add to the reaction (idtype:[value]).
            Defaults to {}.
        - idprefix (str, optional):
            Prefix for the pseudo-identifier. Defaults to 'refineGEMs'.
        - namespace (Literal['BiGG'], optional):
            Namespace to use for the reaction ID.
            If namespace cannot be matched, uses the pseudo-ID
            Defaults to 'BiGG'.

    Returns:
        1. Case successful construction:
                cobra.Reaction:
                    The newly build reaction object.

        2. Case construction not possible:
                None:
                    Nothing to return.

        3. Case reaction found in model.
                list:
                    List of matching reaction IDs (in model).
    """

    # ---------------------
    # check, if ID in model
    # ---------------------
    matches_found = [
        _.id
        for _ in model.reactions
        if "metanetx.reaction" in _.annotation.keys()
        and _.annotation["metanetx.reaction"] == id
    ]
    if len(matches_found) > 0:
        return matches_found

    # -----------------------------
    # otherwise, build new reaction
    # -----------------------------

    # get relevant part of table from database
    mnx_reac_refs = load_a_table_from_database(
        f"SELECT * FROM mnx_reac_xref WHERE id = '{id}'", query=True
    )
    mnx_reac_refs = mnx_reac_refs[
        ~(mnx_reac_refs["description"] == "secondary/obsolete/fantasy identifier")
    ]

    # create reaction object
    new_reac = cobra.Reaction(create_random_id(model, "reac", idprefix))

    # Get names for reaction
    names = []
    reaceqsyms = ["=", "+", "<", ">", "_"]
    important_keywords = ["rna", "dna"]
    for desc in mnx_reac_refs["description"]:
        splitted_descs = re.split(r"\|+", desc)

        # Exclude reaction equations
        names.extend([_ for _ in splitted_descs if not any(s in _ for s in reaceqsyms)])

    # Set name for reaction if possible
    kw_matches = [n for n in names if any(kw in n.lower() for kw in important_keywords)]

    if len(kw_matches) > 0:
        new_reac.name = kw_matches[0]  # Take first entry containing keyword
    elif len(names) > 0:
        new_reac.name = min(names, key=len)  # Take shortest entry
    else:
        new_reac.name = ""

    # get metabolites
    # ---------------
    if not reac_str:
        mnx_reac_prop = load_a_table_from_database(
            f"SELECT * FROM mnx_reac_prop WHERE id = '{id}'", query=True
        )
        reac_str = mnx_reac_prop["mnx_equation"][0]
        if mnx_reac_prop["ec-code"][0]:
            references["ec-code"] = mnx_reac_prop["ec-code"][0]

    if reac_str:
        reactants, products, comparts, rev = parse_reac_str(reac_str, "MetaNetX")
    else:
        return None

    # reac_prop / mnx equation only saves generic compartments 1 and 2 (MNXD1 / MNXD2)
    # current solution for the compartments: 1 -> c, 2 -> e
    comparts = ["c" if _ == "MNXD1" else "e" for _ in comparts]
    metabolites = {}
    meta_counter = 0

    # reconstruct reactants
    for mid, factor in reactants.items():
        tmp_meta = build_metabolite_mnx(
            mid, model, namespace, comparts[meta_counter], idprefix
        )
        if tmp_meta:
            metabolites[tmp_meta] = -1 * factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # reconstruct products
    for mid, factor in products.items():
        tmp_meta = build_metabolite_mnx(
            mid, model, namespace, comparts[meta_counter], idprefix
        )
        if tmp_meta:
            metabolites[tmp_meta] = factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # add metabolites to reaction
    new_reac.add_metabolites(metabolites)

    # set reversibility
    if rev:
        new_reac.bounds = cobra.Configuration().bounds
    else:
        new_reac.bounds = (0.0, cobra.Configuration().upper_bound)

    # get annotations
    # ---------------
    new_reac.annotation["sbo"] = "SBO:0000167"
    # get more annotation from the mnx_reac_xref table
    for db in [
        "bigg.reaction",
        "kegg.reaction",
        "seed.reaction",
        "metacyc.reaction",
        "rhea",
    ]:
        db_matches = mnx_reac_refs[mnx_reac_refs["source"].str.contains(db)]
        if len(db_matches) > 0:
            new_reac.annotation[db] = [
                m.split(":", 1)[1] for m in db_matches["source"].tolist()
            ]
            # update reactions direction, if MetaCyc has better information
            if db == "metacyc.reaction" and len(
                db_matches[db_matches["source"].str.contains("-->")]
            ):
                new_reac.bounds = (0.0, cobra.Configuration().upper_bound)
    # get alias annotations for MetaNetX if available
    if "alias" in references.keys():
        if references["alias"] != None:
            references["metanetx.reaction"] = list(references["alias"])
        del references["alias"]

    # add additional references from the parameter
    _add_annotations_from_dict_cobra(references, new_reac)

    # get annotations
    # ---------------
    new_reac.annotation["sbo"] = "SBO:0000167"

    # add notes
    # ---------
    new_reac.notes["created with"] = "refineGEMs based on MetaNetX"

    # match ID to namespace
    # ---------------------
    match_id_to_namespace(new_reac, namespace)
    # re-check, if reaction not in the model based on the namespace
    if new_reac.id in [_.id for _ in model.reactions]:
        return [new_reac.id]

    return new_reac


def build_reaction_kegg(
    model: cobra.Model,
    id: str = None,
    reac_str: str = None,
    references: dict = {},
    idprefix: str = "refineGEMs",
    namespace: Literal["BiGG"] = "BiGG",
) -> Union[cobra.Reaction, None, list]:
    """Construct a new reaction for a model from either a KEGG reaction ID
    or a KEGG equation string.
    This function will NOT add the reaction directly to the model, if the
    construction process is successful.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - id (str,optional):
            A KEGG reaction ID.
        - reac_str (str, optional):
            The reaction equation string from the database.
            Defaults to None.
        - references (dict, optional):
            Additional annotations to add to the reaction (idtype:[value]).
            Defaults to {}.
        - idprefix (str, optional):
            Prefix for the pseudo-identifier. Defaults to 'refineGEMs'.
        - namespace (Literal['BiGG'], optional):
            Namespace to use for the reaction ID.
            If namespace cannot be matched, uses the pseudo-ID
            Defaults to 'BiGG'.

    Returns:
        1. Case successful construction:
                cobra.Reaction:
                    The newly build reaction object.

        2. Case construction not possible:
                None:
                    Nothing to return.

        3. Case reaction found in model.
                list:
                    List of matching reaction IDs (in model).
    """

    # either reaction id or a reaction string needed for reconstruction
    if not id and not reac_str:
        return None  # reconstruction not possible

    # create an empty reaction with random id
    new_reac = cobra.Reaction(create_random_id(model, "reac", idprefix))

    # -------------
    # KEGG ID given
    # -------------
    if id:
        # check, if reaction in model
        matches = [
            _.id
            for _ in model.reactions
            if "kegg.reaction" in _.annotation.keys()
            and _.annotation["kegg.reaction"] == id
        ]
        if len(matches) > 0:
            return matches  # return matched reaction ids in list

        # retrieve information from KEGG
        kegg_info = kegg_reaction_parser(id)
        if kegg_info:
            if "name" in kegg_info.keys():
                new_reac.name = kegg_info["name"]
            if "equation" in kegg_info.keys():
                reac_str = kegg_info["equation"] if not reac_str else reac_str
            if "db" in kegg_info.keys():
                new_reac.annotation = kegg_info["db"]

    # -------------------------------------
    # Reaction reconstruction from equation
    # -------------------------------------

    # skip, if no reaction string is available
    if not reac_str:
        return None  # reconstruction not possible

    # parse reaction string
    reactants, products, comparts, rev = parse_reac_str(reac_str, "KEGG")

    # KEGG has no information about compartments !!!
    # current solution for the compartments: always use c
    compartment = "c"
    metabolites = {}
    meta_counter = 0

    # reconstruct reactants
    for mid, factor in reactants.items():
        tmp_meta = build_metabolite_kegg(mid, model, namespace, compartment, idprefix)
        if tmp_meta:
            metabolites[tmp_meta] = -1 * factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # reconstruct products
    for mid, factor in products.items():
        tmp_meta = build_metabolite_kegg(mid, model, namespace, compartment, idprefix)
        if tmp_meta:
            metabolites[tmp_meta] = factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # add metabolites to reaction
    new_reac.add_metabolites(metabolites)

    # set reversibility
    if rev:
        new_reac.bounds = (1000.0, 1000.0)
    else:
        new_reac.bounds = (0.0, 1000.0)

    # --------------------
    # add more information
    # --------------------
    new_reac.annotation["sbo"] = "SBO:0000167"
    # get more information from searching the KEGG ID in BiGG
    bigg_res = load_a_table_from_database(
        f"SELECT * FROM bigg_reactions WHERE \"KEGG Reaction\" = '{id}'", query=True
    )
    is_exchange = 0  # check for exchange reactions
    for idx, row in bigg_res.iterrows():
        r, p, compartments, r = parse_reac_str(row["reaction_string"], "BiGG")
        valid_comp = validate_reaction_compartment_bigg(compartments)

        if valid_comp:
            new_reac.annotation["bigg.reaction"] = row["id"]
            _add_annotations_from_bigg_reac_row(row, new_reac)
            if valid_comp == "exchange":
                is_exchange += 1
    new_reac.name = f"EX_{new_reac.name}" if is_exchange > 0 else new_reac.name

    # add additional references from the parameter
    _add_annotations_from_dict_cobra(references, new_reac)

    # add notes
    # ---------
    new_reac.notes["created with"] = "refineGEMs based on KEGG"

    # match ID to namespace
    # ---------------------
    match_id_to_namespace(new_reac, namespace)
    # re-check, if reaction not in the model based on the namespace
    if new_reac.id in [_.id for _ in model.reactions]:
        return [new_reac.id]

    return new_reac


def build_reaction_bigg(
    model: cobra.Model,
    id: str,
    reac_str: str = None,
    references: dict = {},
    idprefix: str = "refineGEMs",
    namespace: Literal["BiGG"] = "BiGG",
) -> Union[cobra.Reaction, None, list]:
    """Construct a new reaction for a model from a BiGG reaction ID.
    This function will NOT add the reaction directly to the model, if the
    construction process is successful.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - id (str):
            A BiGG reaction ID.
        - reac_str (str, optional):
            The reaction equation string from the database.
            Currently, this param is not doing anything in this function.
            Defaults to None.
        - references (dict, optional):
            Additional annotations to add to the reaction (idtype:[value]).
            Defaults to {}.
        - idprefix (str, optional):
            Prefix for the pseudo-identifier. Defaults to 'refineGEMs'.
        - namespace (Literal['BiGG'], optional):
            Namespace to use for the reaction ID.
            If namespace cannot be matched, uses the pseudo-ID
            Defaults to 'BiGG'.

    Returns:
        1. Case successful construction:
                cobra.Reaction:
                    The newly build reaction object.

        2. Case construction not possible:
                None:
                    Nothing to return.

        3. Case reaction found in model.
                list:
                    List of matching reaction IDs (in model).
    """

    # ---------------------
    # check, if ID in model
    # ---------------------
    matches_found = [
        _.id
        for _ in model.reactions
        if "bigg.reaction" in _.annotation.keys()
        and _.annotation["bigg.reaction"] == id
    ]
    if len(matches_found) > 0:
        return matches_found

    # -----------------------------
    # otherwise, build new reaction
    # -----------------------------
    # create reaction object
    new_reac = cobra.Reaction(create_random_id(model, "reac", idprefix))

    # get information from the database
    bigg_reac_info = load_a_table_from_database(
        f"SELECT * FROM bigg_reactions WHERE id = '{id}'", query=True
    ).iloc[0, :]
    new_reac.name = bigg_reac_info["name"]

    # get information from the reaction string
    reactants, products, comparts, rev = parse_reac_str(
        bigg_reac_info["reaction_string"], "BiGG"
    )

    # validate & add compartment information
    # --------------------------------------
    valid_comp = validate_reaction_compartment_bigg(comparts)
    if valid_comp:
        if valid_comp == "exchange":
            new_reac.name = f"EX_{new_reac.name}"
    else:
        return None

    # add metabolites
    # ---------------
    metabolites = {}
    meta_counter = 0
    # reconstruct reactants
    for mid, factor in reactants.items():
        tmp_meta = build_metabolite_bigg(mid, model, namespace, idprefix)
        if tmp_meta:
            metabolites[tmp_meta] = -1 * factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # reconstruct products
    for mid, factor in products.items():
        tmp_meta = build_metabolite_bigg(mid, model, namespace, idprefix)
        if tmp_meta:
            metabolites[tmp_meta] = factor
            meta_counter += 1
        else:
            return None  # not able to build reaction successfully

    # add metabolites to reaction
    new_reac.add_metabolites(metabolites)

    # set reversibility
    if rev:
        new_reac.bounds = (1000.0, 1000.0)
    else:
        new_reac.bounds = (0.0, 1000.0)

    # add annotations
    # ---------------
    # add SBOTerm
    new_reac.annotation["sbo"] = "SBO:0000167"
    # add infos from BiGG
    new_reac.annotation["bigg.reaction"] = [id]
    _add_annotations_from_bigg_reac_row(bigg_reac_info, new_reac)
    # get alias annotations for MetaNetX if available
    if "alias" in references.keys():
        if references["alias"] != None:
            references["metanetx.reaction"] = list(references["alias"])
        del references["alias"]
    # add additional references from the parameter
    _add_annotations_from_dict_cobra(references, new_reac)

    # add notes
    # ---------
    new_reac.notes["created with"] = "refineGEMs based on BiGG"

    # match ID to namespace
    # ---------------------
    match_id_to_namespace(new_reac, namespace)
    # re-check, if reaction not in the model based on the namespace
    if new_reac.id in [_.id for _ in model.reactions]:
        return [new_reac.id]

    return new_reac


@implement
# maybe for later, if we need something to work independantly from namespace
# and outside the gapfilling
# Or would it be better to have a function add_reaction() that only adds
# reactions built with the previuous functions?
def build_reaction():
    pass


# remove genes from model
# -----------------------


def remove_non_essential_genes(
    model: cobra.Model,
    genes_to_check: Union[list[cobra.Gene], None] = None,
    remove_reactions: bool = True,
    min_growth_threshold: float = MIN_GROWTH_THRESHOLD,
) -> tuple[int, int]:
    """Remove non-essential genes based on gene knock-out.

    Args:
        - model (cobra.Model):
            The model to delete the genes from.
        - genes_to_check (Union[list[cobra.Gene],None], optional):
            List of genes to check. If None, all genes of the model are checked.
            Defaults to None.
        - remove_reactions (bool, optional):
            If set to True, also remove corresponding reactions.
            Defaults to True.
        - min_growth_threshold (float, optional):
            Minimal value for growth value.
            Defaults to MIN_GROWTH_THRESHOLD.

    Returns:
        tuple[int,int]:
            Tuple of (1) int, (2) int:
                (1) Number of deleted genes.
                (2) Number of essential genes.
    """

    # initialise values
    to_delete = 0
    essential_counter = 0
    remove = False

    # if no list of genes is given, take all genes of model
    if not genes_to_check:
        genes_to_check = model.genes

    # iterate over all genes of interest
    for g in genes_to_check:
        # test knock out on a model copy
        with model as m:
            # knock out current gene
            g.knock_out()
            res = model.optimize()
            # make sure, model still reaches minimal growth
            if res.objective_value > min_growth_threshold:
                # set for deletion
                remove = True
                to_delete += 1
            else:
                # keep
                essential_counter += 1

        # delete gene
        if remove:
            cobra.manipulation.delete.remove_genes(
                model, [g], remove_reactions=remove_reactions
            )
            remove = False

    return (to_delete, essential_counter)


# ++++++++++++++++++++++++++++++++++++++++
# libSBML - models
# ++++++++++++++++++++++++++++++++++++++++

# extracting reactions & Co via libsbml
# -------------------------------------


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
# Might be possible to deprecate
def compare_gene_lists(
    gps_in_model: pd.DataFrame, db_genes: pd.DataFrame, kegg: bool = True
) -> pd.DataFrame:
    """Compares the provided tables according to column 0/'Locus_tag'

    Args:
        - gps_in_model (pd.DataFrame):
            Table containing the KEGG Gene IDs/Locus tags in the model
        - db_genes (pd.DataFrame):
            Table containing the KEGG Gene IDs for the organism from KEGG/
            locus tags (Accession-2) from BioCyc
        - kegg (bool):
            True if KEGG Genes should be extracted, otherwise False

    Returns:
        pd.DataFrame:
            Table containing all missing genes
    """
    in_db = db_genes.set_index(0) if kegg else db_genes.set_index("locus_tag")
    in_model = gps_in_model.set_index(0)

    genes_in_db_not_in_model = in_db[~in_db.index.isin(in_model.index)]

    return (
        genes_in_db_not_in_model.reset_index().iloc[:, 0]
        if kegg
        else genes_in_db_not_in_model.reset_index()
    )


# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def get_model_reacs_or_metabs(
    model_libsbml: libModel, metabolites: bool = False, col_name: str = "bigg_id"
) -> pd.DataFrame:
    """Extracts table of reactions/metabolites with BiGG IDs from model

    Args:
        - model_libsbml (libModel):
            Model loaded with libSBML
        - metabolites (bool):
            Set to True if metabolites from model should be extracted
        - col_name (str):
            Name to be used for column in table. Defaults to 'bigg_id'.

    Returns:
        pd.DataFrame:
            Table with model identifiers for either metabolites or reactions
    """
    reac_or_metab_list = (
        model_libsbml.getListOfSpecies()
        if metabolites
        else model_libsbml.getListOfReactions()
    )

    list_of_reacs_or_metabs = []
    for reac_or_metab in reac_or_metab_list:
        list_of_reacs_or_metabs.append(reac_or_metab.id[2:])

    reac_or_metab_list_df = pd.Series(list_of_reacs_or_metabs)
    reac_or_metab_list_df = pd.DataFrame(reac_or_metab_list_df, columns=[col_name])

    return reac_or_metab_list_df


def get_reversible(fluxes: dict[str:str]) -> bool:
    """Infer if reaction is reversible from flux bounds

    Args:
        - fluxes (dict):
            Dictionary containing the keys 'lower_bound' & 'upper_bound'
            with values in ['cobra_default_lb', 'cobra_0_bound', 'cobra_default_ub']

    Returns:
        bool:
            True if reversible else False
    """
    return (fluxes["lower_bound"] == "cobra_default_lb") and (
        fluxes["upper_bound"] == "cobra_default_ub"
    )


# map model entities via libSBML
# ------------------------------


def get_gpid_mapping(
    model: libModel,
    gff_paths: Union[str,list[str]] = None,
    email: str = None,
    contains_locus_tags: bool = False,
    outpath: Union[str, Path] = None,
) -> pd.DataFrame:
    """Generate a mapping from model IDs to valid database IDs via model content, GFF files (optional)
    and NCBI requests (optional).

    .. hint:

        Mappings may be incomplete and require manual adjustment.

    .. note:

        Mapping currently is only implemented for NCBI Protein and RefSeq IDs.
        Other databases may be added in the future.

    Args:
        - model (libModel):
            Model loaded with libSBML
        - gff_paths (str|list[str]):
            Path(s) to GFF file(s). Allowed GFF formats are: RefSeq, NCBI and Prokka.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - email (str):
            E-mail for NCBI queries.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - contains_locus_tags (bool, optional):
            Specifies if provided model has locus tags within the label tag if set to True.
            This is only used when mapping_tbl_file == None.
            Defaults to False.
        - outpath (str|Path, optional):
            Output path for location where the generated mapping table should be written to.
            This is only used when mapping_tbl_file == None.
            Defaults to None.

    Returns:
        pd.DataFrame:
            Mapping from model IDs to valid database IDs
    """
    # Set-up path
    if outpath:
        if isinstance(outpath, str):
            outpath = Path(outpath)
        outpath.mkdir(parents=True, exist_ok=True)

    # Helper function to classify potential database IDs according to db specific regexes
    def _classify_potential_db_ids(row: pd.Series) -> pd.Series:
        """Calssifies potential database IDs according to db specific regexes

        Args:
            - row (pd.Series): 
                Row of a DataFrame containing the potential database ID

        Returns:
            pd.Series: 
                Row of the DataFrame with classified database IDs
        """
        # Get NCBI ID
        ncbi_id = row["database_id"] if not pd.isnull(row["database_id"]) else ''

        # Get locus_tag
        locus_tag = row.get("locus_tag", '') if not pd.isnull(row.get("locus_tag")) else ''

        # Check that potential database ID is not the locus_tag
        if (not re.fullmatch(locus_tag, ncbi_id, re.IGNORECASE)) or ('model_id' in row.index.to_list()):
        
            # If identifier matches RefSeq ID pattern
            if re.fullmatch(rf'{DB2REGEX["refseq"]}', ncbi_id, re.IGNORECASE):
                row["REFSEQ"] = row["database_id"]

            # If identifier matches ncbiprotein ID pattern
            elif re.fullmatch(rf'{DB2REGEX["ncbiprotein"]}$', ncbi_id, re.IGNORECASE):
                row["NCBI"] = row["database_id"]

            # No pattern match & no locus_tag match
            else:
                row["UNCLASSIFIED"] = row["database_id"]

        return row

    # 1. Get model IDs & potential contained IDs
    print("Extracting model IDs and potential valid database IDs from model...")
    gene_list = model.getPlugin("fbc").getListOfGeneProducts()
    modelid2potentialid = {"model_id": [], "database_id": []}

    for gene in tqdm(gene_list):

        # Get locus_tag if available
        if contains_locus_tags:
            modelid2potentialid["locus_tag"] = (
                gene.getLabel() if gene.isSetLabel() else None
            )

        # Get current gene id
        gene_id = gene.getId()
        modelid2potentialid["model_id"].append(gene_id)

        # Get gene_id without prefix & replace ascii numbers with according characters
        gene_id = _f_gene(gene_id)

        # Extract potential database ID from gene_id
        if (gene_id.lower() != "spontaneous") and (
            gene_id.lower() != "unknown"
        ):  # Has to be omitted as no additional data can be retrieved neither from NCBI nor the CarveMe input file
            # Case 1: CarveMe version; Contains NCBI Protein ID, lcl_CP035291_1_prot_QCY37541_1_361
            if "prot_" in gene_id:
                 # All NCBI CDS protein FASTA files have the NCBI protein identifier + Version after 'prot_' in the FASTA identifier
                id_string = gene_id.split(r"prot_")[1].split(r"_")
                # If identifier contains no '_', this is full identifier + add version number
                ncbi_id = '.'.join(id_string[0:2])

            # Case 2: RAST contained without any surrounding stuff, 1134914.3.1020_56.peg
            # Case 3: RefSeq contained without any surrounding stuff, WP_011275285.1, From CarveMe: YP_005228568_1
            # Case 4: KEGG contained without any surrounding stuff, sha:SH1836
            # Case 5: Model ID contains locus_tag, AB-1_S128_00983
            else:
                ncbi_id = gene_id

        # Add potential database ID to dictionary
        modelid2potentialid["database_id"].append(ncbi_id)

    # Generate table from dictionary
    mapping_table = pd.DataFrame.from_dict(modelid2potentialid)

    # Add empty column for REFSEQ, NCBI and UNCLASSIFIED IDs
    mapping_table["REFSEQ"] = None
    mapping_table["NCBI"] = None
    mapping_table["UNCLASSIFIED"] = None

    # Classify potential database IDs according to db specific regexes
    print('Classifying potential database IDs according to db specific regexes...')
    mapping_table = mapping_table.progress_apply(_classify_potential_db_ids, axis=1)
    mapping_table.drop("database_id", axis=1, inplace=True)

    # Drop all columns containing only NaNs
    mapping_table = mapping_table.dropna(axis=1, how="all")

    # 1.1 (Optional) Get information from GFF(s)
    if gff_paths:
        print("Extracting (protein id,) locus tag and name from GFF(s)...")

        # Identical attributes to keep per GFF
        to_keep = {"locus_tag": "locus_tag", "product": "name"}

        # Parse & collect all provided GFFs
        gffs = []
        if isinstance(gff_paths, str): gff_paths = [gff_paths]
        for gffp in tqdm(gff_paths):
            current_gff = parse_gff_for_cds(gffp)

            # Check if protein_id in columns
            if "protein_id" in current_gff.columns: to_keep.update({"protein_id": "database_id"})

            # Subset dataframe to only keep relevant columns & rename them
            current_gff = current_gff[to_keep.keys()]
            current_gff = current_gff.rename(columns=to_keep)

            # Explode columns containing lists to get single entries per row
            current_gff = current_gff.explode(column=current_gff.columns.to_list(),ignore_index=True)
            current_gff = current_gff.drop_duplicates()
            current_gff = current_gff.reset_index(drop=True)

            # If protein_id/database_id is available, map it to according database
            if "database_id" in current_gff.columns:

                # Add empty column for REFSEQ, NCBI and UNCLASSIFIED IDs
                current_gff["REFSEQ"] = None
                current_gff["NCBI"] = None
                current_gff["UNCLASSIFIED"] = None

                # Classify potential database IDs according to db specific regexes
                print('Classifying potential database IDs according to db specific regexes...')
                current_gff = current_gff.progress_apply(_classify_potential_db_ids, axis=1)
                current_gff.drop("database_id", axis=1, inplace=True)
                
                # Drop all columns containing only NaNs
                current_gff = current_gff.dropna(axis=1, how="all")

                # Rename columns to avoid duplicates
                if not ('NCBI' in current_gff.columns and 'REFSEQ' in current_gff.columns):
                    if 'NCBI' in current_gff.columns: current_gff.rename(columns={"name": "name_ncbi"}, inplace=True)
                    if 'REFSEQ' in current_gff.columns: current_gff.rename(columns={"name": "name_refseq"}, inplace=True)

            # Add current_gff to GFF list
            gffs.append(current_gff)

        # Merge all GFFs on common column locus_tag
        if len(gffs) == 1:
            gff_mapping = gffs[0]
        else:
            gff_mapping = reduce(
                lambda left, right: pd.merge(
                    left, right, on=["locus_tag"], how="outer", sort=True
                ),
                gffs,
            )

        # Find identical columns in mapping_table & gff_mapping
        identical_columns = list(
            set(mapping_table.columns).intersection(set(gff_mapping.columns))
        )

        # Merge both dataframes
        # If merge on identical columns is possible
        if identical_columns:
            mapping_table = pd.merge(
                mapping_table, gff_mapping, how="outer", on=identical_columns
            )
        else:  
            # Try merge on locus_tag & UNCLASSIFIED as UNCLASSIFIED could contain locus tags
            try:
                mapping_table = pd.merge(
                    mapping_table,
                    gff_mapping,
                    how="outer",
                    left_on="UNCLASSIFIED",
                    right_on="locus_tag",
                )
            except KeyError:
                filename = (
                    f'{model.getId()}_gp_id_gff_mapping_{str(date.today().strftime("%Y%m%d"))}.csv'
                    )
                if outpath:
                    filename = Path(outpath, filename)
                else:
                    filename = Path(filename)
                print(f'''
No common columns found between mapping table derived from model and mapping table derived from provided GFFs. Cannot merge.
The resulting mapping tables will be returned separately. The table for the GFF mapping is written to {filename}.
                ''')
                gff_mapping.to_csv(filename, index=False)

    # 2. Get information from NCBI
    contains_protein_ids = len({'NCBI', 'REFSEQ'}.intersection(set(mapping_table.columns))) > 0
    if email and contains_protein_ids:
        Entrez.email = email
        print("Retrieve locus tag and name information from NCBI...")

        # Add empty columns for names and locus tags if not already contained
        if not "locus_tag" in mapping_table.columns:
            mapping_table["locus_tag"] = None

        contains_only_name = ('name' in mapping_table.columns) and not ('name_refseq' in mapping_table.columns) and not ('name_ncbi' in mapping_table.columns)
        if not any(["name" in _ for _ in mapping_table.columns]) or contains_only_name:
            if ('NCBI' in current_gff.columns and 'REFSEQ' in current_gff.columns):
                mapping_table["name_ncbi"] = None
                mapping_table["name_refseq"] = None
            elif 'NCBI' in current_gff.columns:
                mapping_table["name_ncbi"] = None
            elif 'REFSEQ' in current_gff.columns:
                mapping_table["name_refseq"] = None

        # Query NCBI based on RefSeq Protein ID
        print('Querying NCBI based on RefSeq Protein IDs...')
        if "REFSEQ" in mapping_table.columns:
            mapping_table = mapping_table.progress_apply(
                _search_ncbi_for_gp, axis=1, args=("refseq",)
            )

        # Query NCBI based on NCBI Protein ID
        print('Querying NCBI based on NCBI Protein IDs...')
        if "NCBI" in mapping_table.columns:
            mapping_table = mapping_table.progress_apply(
                _search_ncbi_for_gp, axis=1, args=("ncbiprotein",)
            )

    # Clean-up dataframe for output
    def _merge_name_cols(row: pd.Series) -> pd.Series:
        """Merge all name columns into one column name

        Args:
            - row (pd.Series):
                Row of a DataFrame

        Returns:
            pd.Series:
                DataFrame with only one name column
        """
        # Prefer NCBI Protein name over all other names as most specific for organism
        if 'name_ncbi' in row.index.to_list():
            if row['name_ncbi']:
                row['name'] = row['name_ncbi'] if row['name_ncbi'] != row['NCBI'] else None
            # If no 'valid' NCBI Protein name was found, use REFSEQ name
            elif 'name_refseq' in row.index.to_list():
                if row['name_refseq']:
                    row['name'] = row['name_refseq'] if row['name_refseq'] != row['REFSEQ'] else None
        # If no NCBI Protein name column exists, use REFSEQ name
        elif 'name_refseq' in row.index.to_list():
            if row['name_refseq']:
                row['name'] = row['name_refseq'] if row['name_refseq'] != row['REFSEQ'] else None

        # If no REFSEQ and no NCBI Protein name column exists, check name column for locus_tag
        if row['name'] and (row['name'] == row['locus_tag']):
            row['name'] = None

        # Default name is whatever is stored in name column, if exists, or None
        return row

    print('Cleaning up dataframe for output...')
    # Generate default column name to set content
    if not 'name' in mapping_table.columns:
        mapping_table["name"] = None
    mapping_table = mapping_table.progress_apply(_merge_name_cols, axis=1)
    # Remove unnecessary name columns, if they exist
    names_to_remove = [n_col for n_col in mapping_table.columns if n_col.startswith("name_")]
    if names_to_remove: mapping_table.drop(names_to_remove, axis=1, inplace=True)

    # Write mapping information to file
    filename = (
        f'{model.getId()}_gp_id_mapping_{str(date.today().strftime("%Y%m%d"))}.csv'
    )
    if outpath:
        filename = Path(outpath, filename)
    else:
        filename = Path(filename)
    logging.info(
        f"The mapping table mapping the geneProduct IDs of the provided model {model.getId()} to the names and locus "
        + f"tags is saved to {filename}"
    )
    mapping_table.to_csv(filename, index=False)

    return mapping_table


# create model entities using libSBML
# -----------------------------------

def create_gp(
    model: libModel,
    protein_id: str,
    model_id: str = None,
    name: str = None,
    locus_tag: str = None,
    reference: dict[str : tuple[Union[list, str], bool]] = dict(),
    sanity_check: bool = False
) -> None:
    """Creates GeneProduct in the given libSBML model.

    Args:
        - model (libModel):
            The model object, loaded with libSBML.
        - protein_id (str):
            (NCBI) Protein ID of the gene.
        - model_id (str, optional):
            If given, uses this string as the ID of the gene in the model.
            ID should be identical to ID that CarveMe adds from the NCBI FASTA input file.
            Defaults to None.
        - name (str, optional):
            Name of the GeneProduct. Defaults to None.
        - locus_tag (str, optional):
            Genome-specific locus tag. Will be used as label in the model.
            Defaults to None.
        - reference (dict, optional):
            Dictionary containing references for the gene product.
            The key is the database name, the value is a tuple with the first element being the ID(s)
            and the second element being a boolean indicating if the strain is a lab strain or not.
            Defaults to an empty dictionary.
        - sanity_check (bool, optional):
            Check, whether locus tag (label) or model ID (ID) already exist in model.
            Note, that setting this to True increases the runtime.
            Defaults to False.
    """

    # check, if possible ID has been passed
    if model_id:  # ID
       pass
    elif protein_id:
        model_id = _f_gene_rev(protein_id)
    elif locus_tag:
        model_id = _f_gene_rev(locus_tag)
    else:
        raise ValueError(
            "No valid ID found for gene product. Specify at least the locus tag."
        )
    
    # construct label, if possible
    label = None
    if locus_tag:
        label = _f_gene_rev(locus_tag, "")
    
    # sanity check
    if sanity_check:
        genes = model.getPlugin(0).getListOfGeneProducts()
        for g in genes:
            # for ID 
            if g.isSetId() and g.getId == model_id:
                logging.warning(f'Cannot add gene product, as ID {model_id} is alreaedy in the model.')
                return None
            # for locus tag as label
            if label and g.isSetLabel() and g.getLabel() == label:
                logging.warning(f'Cannot add gene product, as label {label} is already in the model.')
                return None
        
    # create gene product object
    gp = model.getPlugin(0).createGeneProduct()
    # set basic attributes
    gp.setIdAttribute(model_id)
    if name:
        gp.setName(name)  # Name
    if label:
        gp.setLabel(label)  # Label
    gp.setSBOTerm("SBO:0000243")  # SBOterm
    gp.setMetaId(f"meta_{_f_gene_rev(protein_id)}")  # Meta ID
    # test for NCBI/RefSeq
    id_db = None
    if re.fullmatch(
        r"^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$",
        protein_id,
        re.IGNORECASE,
    ):
        id_db = "REFSEQ"
    elif re.fullmatch(r"^(\w+\d+(\.\d+)?)|(NP_\d+)$", protein_id, re.IGNORECASE):
        id_db = "NCBI"
    if id_db:
        add_cv_term_genes(protein_id, id_db, gp)  # NCBI protein

    # add further references
    # references {'dbname':(list_or_str_of_id,lab_strain_bool)}
    if reference and len(reference) > 0:
        for dbname, specs in reference.items():
            match specs[0]:
                case list():
                    for s in specs[0]:
                        add_cv_term_genes(s, REF_COL_GF_GENE_MAP[dbname], gp, specs[1])
                case str():
                    add_cv_term_genes(
                        specs[0], REF_COL_GF_GENE_MAP[dbname], gp, specs[1]
                    )
                case None:
                    pass
                case _:
                    raise ValueError(f"Unexpected type for reference value: {specs[0]}")


def create_species(
    model: libModel,
    metabolite_id: str,
    name: str,
    compartment_id: str,
    charge: int,
    chem_formula: str,
) -> tuple[Species, libModel]:
    """Creates Species/Metabolite in the given model

    Args:
        - model (libModel):
            Model loaded with libSBML
        - metabolite_id (str):
            Metabolite ID within model (If model from CarveMe, preferable a BiGG ID)
        - name (str):
            Name of the metabolite
        - compartment_id (str):
            ID of the compartment where metabolite resides
        - charge (int):
            Charge for the metabolite
        - chem_formula (str):
            Chemical formula for the metabolite

    Returns:
        tuple:
            libSBML Species (1) & libSBML model (2)

            (1) Species: Created species/metabolite
            (2) libModel: Model containing the created metabolite
    """
    metabolite = model.createSpecies()
    metabolite.setId(f"M_{metabolite_id}")
    if name:
        metabolite.setName(name)
    metabolite.setMetaId(f"meta_M_{metabolite_id}")
    metabolite.setSBOTerm("SBO:0000247")
    metabolite.setInitialAmount(float("NaN"))
    metabolite.setHasOnlySubstanceUnits(True)
    metabolite.setBoundaryCondition(False)
    metabolite.setConstant(False)
    metabolite.setCompartment(compartment_id)
    metabolite.getPlugin(0).setCharge(charge)
    metabolite.getPlugin(0).setChemicalFormula(chem_formula)
    add_cv_term_metabolites(metabolite_id[:-2], "BIGG", metabolite)
    return metabolite, model


def create_reaction(
    model: libModel,
    reaction_id: str,
    name: str,
    reactants: dict[str:int],
    products: dict[str:int],
    fluxes: dict[str:str],
    reversible: bool = None,
    fast: bool = None,
    compartment: str = None,
    sbo: str = None,
    genes: Union[str, list[str]] = None,
) -> tuple[Reaction, libModel]:
    """Creates new reaction in the given model

    Args:
        - model (libModel):
            Model loaded with libSBML
        - reaction_id (str):
            BiGG ID of the reaction to create
        - name (str):
            Human readable name of the reaction
        - reactants (dict):
            Metabolites as keys and their stoichiometry as values
        - products (dict):
            Metabolites as keys and their stoichiometry as values
        - fluxes (dict):
            Dictionary with lower_bound and upper_bound as keys
        - reversible (bool):
            True/False for the reaction
        - fast (bool):
            True/False for the reaction
        - compartment (str):
            BiGG compartment ID of the reaction (if available)
        - sbo (str):
            SBO term of the reaction
        - genes (str|list):
            List of genes belonging to reaction

    Returns:
        tuple:
            libSBML reaction (1) & libSBML model (2)

            (1) Reaction: Created reaction
            (2) libModel: Model containing the created reaction
    """
    reaction = model.createReaction()
    reaction.setId("R_" + reaction_id)
    if name:
        reaction.setName(name)
    reaction.setMetaId("meta_R_" + reaction_id)
    sbo = (
        sbo if sbo else "SBO:0000167"
    )  # SBO term for biochemical or transport reaction
    reaction.setSBOTerm(sbo)
    fast = fast if fast else False
    reaction.setFast(fast)
    if compartment:
        reaction.setCompartment(
            compartment
        )  # Set compartment for reaction if available
    reversible = reversible if reversible else get_reversible(fluxes)
    reaction.setReversible(reversible)
    if genes:
        if genes == "G_spontaneous":
            reaction.getPlugin(
                0
            ).createGeneProductAssociation().createGeneProductRef().setGeneProduct(gene)
        elif len(genes) == 1:
            reaction.getPlugin(
                0
            ).createGeneProductAssociation().createGeneProductRef().setGeneProduct(
                genes[0]
            )
        else:
            gp_ass_or = reaction.getPlugin(0).createGeneProductAssociation().createOr()
            for gene in genes:
                # Set GeneProductReferences if available
                gp_ass_or.createGeneProductRef().setGeneProduct(gene)
    for metab, stoich in reactants.items():  # reactants as dict with metab:stoich
        reaction.addReactant(model.getSpecies("M_" + metab), stoich)
    for metab, stoich in products.items():  # reactants as dict with metab:stoich
        reaction.addProduct(model.getSpecies("M_" + metab), stoich)
    reaction.getPlugin(0).setLowerFluxBound(fluxes["lower_bound"])
    reaction.getPlugin(0).setUpperFluxBound(fluxes["upper_bound"])
    add_cv_term_reactions(reaction_id, "BIGG", reaction)
    return reaction, model


def create_gpr(reaction: Reaction, gene: Union[str, list[str]]) -> None:
    """For a given libSBML Reaction and a gene ID or a list of gene IDs,
    create a gene production rule inside the reaction.

    Currently only supports 'OR' causality.

    Args:
        - reaction (libsbml.Reaction):
            The reaction object to add the GPR to.
        - gene (str | list[str]):
            Either a gene ID or a list of gene IDs, that will be added to the GPR
            (OR causality).
    """

    # Step 1: test, if there is already a gpr
    # ---------------------------------------
    old_association_str = None
    old_association_fbc = None
    if reaction.getPlugin(0).getGeneProductAssociation():
        old_association = (
            reaction.getPlugin(0).getGeneProductAssociation().getListOfAllElements()
        )
        # case 1: only a single association
        if len(old_association) == 1 and isinstance(old_association[0], GeneProductRef):
            old_association_str = old_association[0].getGeneProduct()
        # case 2: nested structure of asociations
        elif isinstance(old_association[0], FbcOr) or isinstance(
            old_association[0], FbcAnd
        ):
            old_association_fbc = old_association[0].clone()
            # this should get the highest level association (that includes all others)

    # Step 2: create new gene product association
    # -------------------------------------------
    if old_association_str and isinstance(gene, str):
        gene = [old_association_str, gene]
    elif old_association_str and isinstance(gene, list):
        gene.append(old_association_str)

    # add the old association rule as an 'OR' (if needed)
    if not old_association_fbc:
        new_association = reaction.getPlugin(0).createGeneProductAssociation()
    else:
        new_association = (
            reaction.getPlugin(0).createGeneProductAssociation().createOr()
        )
        new_association.addAssociation(old_association_fbc)

    # add the remaining genes
    if isinstance(gene, str):
        new_association.createGeneProductRef().setGeneProduct(gene)
    elif isinstance(gene, list) and len(gene) == 1:
        new_association.createGeneProductRef().setGeneProduct(gene[0])
    elif isinstance(gene, list) and len(gene) > 1:
        gpa_or = new_association.createOr()
        for i in gene:
            gpa_or.createGeneProductRef().setGeneProduct(i)


def create_unit(
    model_specs: tuple[int],
    meta_id: str,
    kind: str,
    e: int,
    m: int,
    s: int,
    uri_is: str = "",
    uri_idf: str = "",
) -> Unit:
    """Creates unit for SBML model according to arguments

    Args:
        - model_specs (tuple):
            Level & Version of SBML model
        - meta_id (str):
            Meta ID for unit (Neccessary for URI)
        - kind (str):
            Unit kind constant (see libSBML for available constants)
        - e (int):
            Exponent of unit
        - m (int):
            Multiplier of unit
        - s (int):
            Scale of unit
        - uri_is (str):
            URI supporting the specified unit
        - uri_idf (str):
            URI supporting the derived from unit

    Returns:
        Unit:
            libSBML unit object
    """
    unit = Unit(*model_specs)
    unit.setKind(kind)
    unit.setMetaId(f"meta_{meta_id}_{unit.getKind()}")
    unit.setExponent(e)
    unit.setMultiplier(m)
    unit.setScale(s)
    if uri_is:
        add_cv_term_units(uri_is, unit, BQM_IS)
    if uri_idf:
        add_cv_term_units(uri_idf, unit, BQM_IS_DERIVED_FROM)
    return unit


def create_unit_definition(
    model_specs: tuple[int], identifier: str, name: str, units: list[Unit]
) -> UnitDefinition:
    """Creates unit definition for SBML model according to arguments

    Args:
        - model_specs (tuple):
            Level & Version of SBML model
        - identifier (str):
            Identifier for the defined unit
        - name (str):
            Full name of the defined unit
        - units (list):
            All units the defined unit consists of

    Returns:
        UnitDefinition:
            libSBML unit definition object
    """
    unit_definition = UnitDefinition(*model_specs)
    unit_definition.setId(identifier)
    unit_definition.setMetaId(f"meta_{identifier}")
    unit_definition.setName(name)

    # Iterate over all units provided for the unit definition & add the units
    for unit in units:
        unit_definition.addUnit(unit)

    return unit_definition


# create default fba units
# ^^^^^^^^^^^^^^^^^^^^^^^^


def create_fba_units(model: libModel) -> list[UnitDefinition]:
    """Creates all fba units required for a constraint-based model

    Args:
        - model (libModel):
            Model loaded with libSBML

    Returns:
        list:
            List of libSBML UnitDefinitions
    """
    # Get model level & version for unit & unit definition
    model_specs = model.getLevel(), model.getVersion()

    # Create required units
    litre = create_unit(
        model_specs,
        "litre",
        UNIT_KIND_LITRE,
        e=1,
        m=1,
        s=-3,
        uri_is="0000104",
        uri_idf="0000099",
    )
    mole = create_unit(
        model_specs,
        "mole_0",
        UNIT_KIND_MOLE,
        e=1,
        m=1,
        s=-3,
        uri_is="0000040",
        uri_idf="0000013",
    )
    gram = create_unit(
        model_specs, "per_gram_0", UNIT_KIND_GRAM, e=-1, m=1, s=0, uri_is="0000021"
    )
    second = create_unit(
        model_specs,
        "second_0",
        UNIT_KIND_SECOND,
        e=1,
        m=3600,
        s=0,
        uri_is="0000032",
        uri_idf="0000010",
    )

    # Create unit definitions for hour & femto litre
    hour = create_unit_definition(
        model_specs, identifier="h", name="Hour", units=[second]
    )
    femto_litre = create_unit_definition(
        model_specs, identifier="fL", name="Femto litres", units=[litre]
    )

    # Create unit definitions for millimoles per gram dry weight (mmgdw) & mmgdw per hour
    mmgdw = create_unit_definition(
        model_specs,
        identifier="mmol_per_gDW",
        name="Millimoles per gram (dry weight)",
        units=[mole, gram],
    )

    # Create new units mole & gram to get new meta IDs
    mole = create_unit(
        model_specs,
        "mole_1",
        UNIT_KIND_MOLE,
        e=1,
        m=1,
        s=-3,
        uri_is="0000040",
        uri_idf="0000013",
    )
    gram = create_unit(
        model_specs, "per_gram_1", UNIT_KIND_GRAM, e=-1, m=1, s=0, uri_is="0000021"
    )
    # Create new unit second to fit to per hour & get new meta ID
    second = create_unit(
        model_specs,
        "second_1",
        UNIT_KIND_SECOND,
        e=-1,
        m=3600,
        s=0,
        uri_is="0000032",
        uri_idf="0000010",
    )

    mmgdwh = create_unit_definition(
        model_specs,
        identifier="mmol_per_gDW_per_h",
        name="Millimoles per gram (dry weight) per hour",
        units=[mole, gram, second],
    )
    add_cv_term_units("pubmed:7986045", mmgdwh, BQM_IS_DESCRIBED_BY)

    return [mmgdwh, mmgdw, hour, femto_litre]


# print model entities using libsbml
# ----------------------------------


def print_UnitDefinitions(unit_defs: ListOfUnitDefinitions):
    """Prints a list of libSBML UnitDefinitions as XMLNodes

    Args:
        - unit_defs (ListOfUnitDefinitions):
            List of libSBML UnitDefinition objects
    """
    for unit_def in unit_defs:
        logging.info(unit_def.toXMLNode())
