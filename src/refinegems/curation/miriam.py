#!/usr/bin/env python
"""General functions to conform annotations to the MIRIAM standards

The functions can be used to change the CURIE pattern and to clean the CURIEs, CVTerm qualifiers and qualifier types
up.
"""

__author__ = "Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import re
import libsbml
import logging

from bioregistry import (
    is_valid_curie,
    is_valid_identifier,
    manager,
    normalize_parsed_curie,
    get_identifiers_org_iri,
)
from colorama import Fore, Style
from datetime import date
from functools import reduce

from libsbml import Model as libModel
from libsbml import UnitDefinition, SBase
from libsbml import (
    MODEL_QUALIFIER,
    BQM_IS,
    BIOLOGICAL_QUALIFIER,
    BQB_IS,
    BQB_HAS_PROPERTY,
    BQB_IS_HOMOLOG_TO,
)
from libsbml import BiolQualifierType_toString, ModelQualifierType_toString

from pathlib import Path
from sortedcontainers import SortedDict, SortedSet
from tqdm.auto import tqdm


from ..utility.cvterms import generate_cvterm, MIRIAM, OLD_MIRIAM
from ..utility.db_access import BIOCYC_TIER1_DATABASES_PREFIXES
from ..utility.io import parse_dict_to_dataframe

################################################################################
# functions
################################################################################


# Change CURIE pattern/Correct CURIEs
# ------------------------------------
def get_set_of_curies(
    uri_list: list[str],
) -> tuple[SortedDict[str : SortedSet[str]], list[str]]:
    """Gets a list of URIs & maps the database prefixes to their respective identifier sets

    Args:
        - uri_list (list[str]):
            List containing CURIEs

    Returns:
        tuple:
            A sorted dictionary (1) & a list (2)

            (1) SortedDict: Sorted dictionary mapping database prefixes from the provided CURIEs to their respective identifier sets also provided by the CURIEs
            (2) list: List of CURIEs that are invalid according to bioregistry
    """
    curie_dict = SortedDict()
    prefix, identifier = None, None
    invalid_curies = []

    for uri in uri_list:

        # Extracts the CURIE part from the URI/IRI
        if MIRIAM in uri:
            extracted_curie = uri.split(MIRIAM)[1]
        elif OLD_MIRIAM in uri:
            extracted_curie = uri.split(OLD_MIRIAM)[1]
        else:
            extracted_curie = uri

        curie = manager.parse_curie(
            extracted_curie
        )  # Contains valid db prefix to identifiers pairs
        curie = list(curie)  # Turn tuple into list to allow item assignment

        # Prefix is valid but to have same result for same databases need to do a bit of own parsing
        if curie[0]:  
            if re.fullmatch(
                r"^biocyc$", curie[0], re.IGNORECASE
            ):  # Check for biocyc to also add metacyc if possible
                # Always add META if BioCyc sub-datbase prefixes are missing
                curie = (
                    curie
                    if curie[1].split(r":")[0] in BIOCYC_TIER1_DATABASES_PREFIXES
                    else [curie[0], f"META:{curie[1]}"]
                )

                if "META" in curie[1]:
                    if is_valid_identifier(
                        *curie
                    ):  # Get the valid BioCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f"{curie[0]}:{curie[1]}")

                    # Add the MetaCyc identifier additionally
                    curie[1] = curie[1].split(r"META:")[
                        1
                    ]  # Metacyc identifier comes after 'META:' in biocyc identifier
                    if re.search(r"^rxn-|-rxn$", curie[1], re.IGNORECASE):
                        curie[0] = "metacyc.reaction"
                    else:
                        curie[0] = "metacyc.compound"
            elif "metacyc." in curie[0]:
                if is_valid_identifier(
                    *curie
                ):  # Get the valid MetaCyc identifier & Add to dictionary
                    prefix, identifier = normalize_parsed_curie(*curie)

                    if not curie_dict or (prefix not in curie_dict):
                        curie_dict[prefix] = SortedSet()
                    curie_dict[prefix].add(identifier)
                else:
                    invalid_curies.append(f"{curie[0]}:{curie[1]}")

                # Add the BioCyc identifier additionally
                curie = [
                    "biocyc",
                    f"META:{curie[1]}",
                ]  # Metacyc identifier comes after 'META:' in biocyc identifier
            elif re.fullmatch(
                r"^brenda$", curie[0], re.IGNORECASE
            ):  # Brenda & EC code is the same
                curie[0] = "eccode"

        elif not curie[0]:  # Need to do own parsing if prefix is not valid
            # Get CURIEs irrespective of pattern
            if "/" in extracted_curie:
                extracted_curie = extracted_curie.split(r"/")

                # Check for NaN identifiers
                if re.fullmatch(
                    r"^nan$", extracted_curie[0], re.IGNORECASE
                ) or re.fullmatch(r"^nan$", extracted_curie[1], re.IGNORECASE):
                    # Only return strings where the database prefix is 'NaN' but a possible identifier could be contained
                    if re.fullmatch(
                        r"^nan$", extracted_curie[0], re.IGNORECASE
                    ) and not re.fullmatch(r"^nan$", extracted_curie[1], re.IGNORECASE):
                        invalid_curies.append(
                            f"{extracted_curie[0]}:{extracted_curie[1]}"
                        )
                    continue
                # Check for certain special cases
                if re.search(
                    r"inchi", extracted_curie[0], re.IGNORECASE
                ):  # Check for inchi as splitting by '/' splits too much
                    if re.fullmatch(r"^inchi$", extracted_curie[0], re.IGNORECASE):
                        curie = (
                            extracted_curie[0].lower(),
                            "/".join(extracted_curie[1 : len(extracted_curie)]),
                        )
                    elif re.fullmatch(r"^inchikey$", extracted_curie[0], re.IGNORECASE):
                        curie = (extracted_curie[0].lower(), extracted_curie[1])
                    else:
                        wrong_prefix = extracted_curie[0].split(r":")
                        curie = (
                            wrong_prefix[0],
                            f'{wrong_prefix[1]}/{"/".join(extracted_curie[1:len(extracted_curie)])}',
                        )
                elif re.fullmatch(
                    r"^brenda$", extracted_curie[0], re.IGNORECASE
                ):  # Brenda & EC code is the same
                    curie = ("eccode", extracted_curie[1])
                elif re.fullmatch(
                    r"^biocyc$", extracted_curie[0], re.IGNORECASE
                ):  # Check for biocyc to also add metacyc if possible
                    # Always add META if BioCyc sub-datbase prefixes are missing
                    extracted_curie[1] = (
                        extracted_curie[1]
                        if extracted_curie[1].split(r":")[0]
                        in BIOCYC_TIER1_DATABASES_PREFIXES
                        else f"META:{extracted_curie[1]}"
                    )
                    curie = ["biocyc", extracted_curie[1]]

                    if "META" in curie[1]:
                        if is_valid_identifier(
                            *curie
                        ):  # Get the valid BioCyc identifier & Add to dictionary
                            prefix, identifier = normalize_parsed_curie(*curie)

                            if not curie_dict or (prefix not in curie_dict):
                                curie_dict[prefix] = SortedSet()
                            curie_dict[prefix].add(identifier)
                        else:
                            invalid_curies.append(f"{curie[0]}:{curie[1]}")

                        # Add additionallly the MetaCyc identifier
                        curie[1] = curie[1].split(r"META:")[
                            1
                        ]  # Metacyc identifier comes after 'META:' in biocyc identifier
                        if re.search(r"^rxn-|-rxn$", curie[1], re.IGNORECASE):
                            curie[0] = "metacyc.reaction"
                        else:
                            curie[0] = "metacyc.compound"
                elif "metacyc." in extracted_curie[0]:
                    curie = extracted_curie
                    if is_valid_identifier(
                        *curie
                    ):  # Get the valid MetaCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f"{curie[0]}:{curie[1]}")

                    # Add BioCyc identfier additionally
                    curie = [
                        "biocyc",
                        f"META:{curie[1]}",
                    ]  # Metacyc identifier comes after 'META:' in biocyc identifier
                elif re.fullmatch(r"^eco|chebi$", extracted_curie[0], re.IGNORECASE):
                    if ":" in extracted_curie[1]:
                        new_curie = extracted_curie[1].split(r":")
                        curie = (new_curie[0].lower(), new_curie[1])
                    else:
                        curie = tuple(extracted_curie)
                elif re.search(
                    r"^sbo:", extracted_curie[1], re.IGNORECASE
                ):  # Checks for old pattern of SBO term URIs ('MIRIAM/sbo/SBO:identifier')
                    curie = [extracted_curie[0], extracted_curie[1].split(r":")[1]]
                else:
                    if re.fullmatch(
                        r"^brenda$", extracted_curie[0], re.IGNORECASE
                    ) or re.fullmatch(
                        r"^ec-code$", extracted_curie[0], re.IGNORECASE
                    ):  # Brenda equals EC code, EC code in URI = ec-code
                        curie[0] = "eccode"
                    else:
                        curie[0] = extracted_curie[0]

                    curie[1] = extracted_curie[1]

            elif ":" in extracted_curie:
                extracted_curie = extracted_curie.split(r":")

                # Check for NaN identifiers
                if re.fullmatch(
                    r"^nan$", extracted_curie[0], re.IGNORECASE
                ) or re.fullmatch(r"^nan$", extracted_curie[1], re.IGNORECASE):
                    # Only return strings where the database prefix is 'NaN' but a possible identifier could be contained
                    if re.fullmatch(
                        r"^nan$", extracted_curie[0], re.IGNORECASE
                    ) and not re.fullmatch(r"^nan$", extracted_curie[1], re.IGNORECASE):
                        invalid_curies.append(
                            f"{extracted_curie[0]}:{extracted_curie[1]}"
                        )
                    continue
                elif re.fullmatch(
                    r"^biocyc$", extracted_curie[0], re.IGNORECASE
                ):  # Check for biocyc to also add metacyc if possible
                    # Always add META if BioCyc sub-datbase prefixes are missing
                    extracted_curie[1] = (
                        extracted_curie[1]
                        if extracted_curie[1].split(r":")[0]
                        in BIOCYC_TIER1_DATABASES_PREFIXES
                        else f"META:{extracted_curie[1]}"
                    )
                    curie = ["biocyc", extracted_curie[1]]

                    if "META" in curie[1]:
                        if is_valid_identifier(
                            *curie
                        ):  # Get the valid BioCyc identifier & Add to dictionary
                            prefix, identifier = normalize_parsed_curie(*curie)

                            if not curie_dict or (prefix not in curie_dict):
                                curie_dict[prefix] = SortedSet()
                            curie_dict[prefix].add(identifier)
                        else:
                            invalid_curies.append(f"{curie[0]}:{curie[1]}")

                        # Add MetaCyc identifier additionally
                        curie[1] = curie[1].split(r"META:")[
                            1
                        ]  # Metacyc identifier comes after 'META:' in biocyc identifier
                        if re.search(r"^rxn-|-rxn$", curie[1], re.IGNORECASE):
                            curie[0] = "metacyc.reaction"
                        else:
                            curie[0] = "metacyc.compound"
                elif "metacyc." in extracted_curie[0]:
                    curie = extracted_curie
                    if is_valid_identifier(
                        *curie
                    ):  # Get the valid MetaCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f"{curie[0]}:{curie[1]}")

                    # Add BioCyc identifier additionally
                    curie = [
                        "biocyc",
                        f"META:{curie[1]}",
                    ]  # Metacyc identifier comes after 'META:' in biocyc identifier
                else:
                    if re.fullmatch(
                        r"^brenda$", extracted_curie[0], re.IGNORECASE
                    ) or re.fullmatch(
                        r"^ec-code$", extracted_curie[0], re.IGNORECASE
                    ):  # Brenda equals EC code, EC code in URI = ec-code
                        curie[0] = "eccode"
                    else:
                        curie[0] = extracted_curie[0]

                    if re.fullmatch(r"^kegg.genes$", extracted_curie[0], re.IGNORECASE):
                        curie[1] = ":".join(extracted_curie[1 : len(extracted_curie)])
                    else:
                        curie[1] = extracted_curie[1]

        if is_valid_identifier(*curie):  # Get all valid identifiers
            prefix, identifier = normalize_parsed_curie(*curie)
        else:

            if curie[0] == "eccode":
                correct_id = curie[
                    1
                ]  # EC number needs to have 4 places if splitted at the dots
                while len(correct_id.split(r".")) < 4:
                    correct_id = f"{correct_id}.-"
                prefix, identifier = normalize_parsed_curie(curie[0], correct_id)
                # Report too long EC codes as invalid CURIEs!
                if len(correct_id.split(r".")) > 4:
                    invalid_curies.append(f"{prefix}:{identifier}")
                    continue
            else:
                invalid_curies.append(f"{curie[0]}:{curie[1]}")

        if prefix and identifier:  # Check that a prefix & identifier pair was found!
            # Use prefix as key & the corresponding set of identifiers as values
            if not curie_dict or (prefix not in curie_dict):
                curie_dict[prefix] = SortedSet()
            curie_dict[prefix].add(identifier)

    return curie_dict, invalid_curies


def generate_uri_set_with_old_pattern(
    prefix2id: SortedDict[str : SortedSet[str]],
) -> SortedSet[str]:
    """Generate a set of complete URIs from the provided prefix to identifier mapping with the old
    MIRIAM pattern.

    Args:
        - prefix2id (SortedDict[str: SortedSet[str]]):
            Dictionary containing a mapping from database prefixes to their respective identifier sets

    Returns:
        SortedSet:
            Sorted set containing complete URIs
    """
    uri_set = SortedSet()

    SEPARATOR = "/"

    for prefix in prefix2id:
        current_prefix = prefix

        for identifier in prefix2id.get(prefix):
            separator = SEPARATOR

            if re.search(
                r"o$", prefix, re.IGNORECASE
            ):  # Ontologies seem only to work with new pattern!
                separator = ":"
                prefix = prefix.upper()

            elif re.fullmatch(
                r"^chebi$", current_prefix, re.IGNORECASE
            ):  # The old pattern for chebi is different: Just adding '/' does NOT work!
                prefix = f"chebi/{current_prefix}"
                separator = ":"

            elif re.fullmatch(
                r"^biocyc$", prefix, re.IGNORECASE
            ):  # Get identifier for biocyc
                prefix = f"biocyc"  # META
                # separator = ':'

            uri = MIRIAM + prefix + separator + identifier
            uri_set.add(uri)

    return uri_set


def generate_miriam_compliant_uri_set(
    prefix2id: SortedDict[str : SortedSet[str]],
) -> SortedSet[str]:
    """Generate a set of complete MIRIAM compliant URIs from the provided prefix to identifier mapping

    Args:
        - prefix2id (SortedDict[str: SortedSet[str]]):
            Dictionary containing a mapping from database prefixes to their respective identifier sets

    Returns:
        SortedSet:
            Sorted set containing complete URIs
    """
    uri_set = SortedSet()

    for prefix in prefix2id:
        for identifier in prefix2id.get(prefix):
            uri = get_identifiers_org_iri(prefix, identifier)
            uri_set.add(uri)

    return uri_set


def add_uri_set(entity: SBase, qt, b_m_qt, uri_set: SortedSet[str]) -> list[str]:
    """Add a complete URI set to the provided CVTerm

    Args:
        - entity (SBase):
            A libSBML SBase object like model, GeneProduct, etc.
        - qt:
            A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - b_m_qt:
            A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - uri_set (SortedSet[str]):
            SortedSet containing URIs
    """
    new_cvterm = generate_cvterm(qt, b_m_qt)

    for uri in uri_set:
        new_cvterm.addResource(uri)

    entity.addCVTerm(new_cvterm)


def improve_uri_per_entity(entity: SBase, new_pattern: bool) -> list[str]:
    """Helper function: Removes duplicates & changes pattern according to new_pattern

    Args:
        - entity (SBase):
            A libSBML SBase object, either a model or an entity
        - new_pattern (bool):
            True if new pattern is wanted, otherwise False

    Returns:
        list:
            List of all collected invalid CURIEs of one entity
    """
    collected_invalid_curies = []
    pattern = rf"{MIRIAM}|{OLD_MIRIAM}"
    cvterms = entity.getCVTerms()
    cvterms = [cvterm.clone() for cvterm in cvterms]
    entity.unsetCVTerms()

    for cvterm in cvterms:
        tmp_list = []
        current_b_m_qt = None  # Needs to be initialised, otherwise UnboundLocalError: local variable 'current_b_m_qt' referenced before assignment

        # Retrieve QualifierType & Biological/ModelQualifierType before resource is removed!
        current_qt = cvterm.getQualifierType()

        match current_qt:

            case libsbml.BIOLOGICAL_QUALIFIER:
                current_b_m_qt = cvterm.getBiologicalQualifierType()

            case libsbml.MODEL_QUALIFIER:
                current_b_m_qt = cvterm.getModelQualifierType()

            case _:
                mes = f"Unknown qualifier type detected in model: {current_qt}"
                raise TypeError(mes)

        current_uris = [
            cvterm.getResourceURI(i) for i in range(cvterm.getNumResources())
        ]

        # Remove all URIs/CURIEs from model & collect invalid/not miriam compliant URIs
        for cu in current_uris:
            current_uri = cu
            cvterm.removeResource(cu)

            if (
                is_valid_identifier(*manager.parse_uri(current_uri))
                or is_valid_curie(current_uri.lower())
                or re.match(pattern, current_uri, re.IGNORECASE)
            ):
                # For Rhea entries if version is specified with '#' remove the version
                if re.search(r"rhea", current_uri, re.IGNORECASE) and re.search(
                    r"#", current_uri
                ):
                    current_uri = current_uri.split(r"#")[0]

                # Collect all URIs to be adjusted/newly added
                tmp_list.append(current_uri)

            else:
                collected_invalid_curies.append(current_uri)

        prefix2id, invalid_curies = get_set_of_curies(tmp_list)
        collected_invalid_curies.extend(invalid_curies)
        if prefix2id:
            if new_pattern:
                uri_set = generate_miriam_compliant_uri_set(prefix2id)
            else:
                uri_set = generate_uri_set_with_old_pattern(prefix2id)
            add_uri_set(entity, current_qt, current_b_m_qt, uri_set)
        else:
            # Remove annotations if no valid URIs/CURIEs were found
            if cvterm.getNumResources() < 1:
                logging.warning(
                    f"No valid URIs/CURIEs found for {entity.getId()}. To resolve manually please inspect file containing invalid CURIEs."
                )

    return collected_invalid_curies


def improve_uris(entities: SBase, new_pattern: bool) -> dict[str : list[str]]:
    """Removes duplicates & changes pattern according to bioregistry or new_pattern

    Args:
        - entities (SBase):
            A libSBML SBase object, either a model or a list of entities
        - bioregistry (bool):
            Specifies whether the URIs should be changed with the help of bioregistry to be MIRIAM compliant or changed according to new or old pattern
        - new_pattern (bool):
            True if new pattern is wanted, otherwise False

    Returns:
        dict:
            Mapping of entity identifier to list of corresponding invalid CURIEs
    """
    entity2invalid_curies = {}

    if type(entities) == libModel:  # Model needs to be handled like entity!
        invalid_curies = improve_uri_per_entity(entities, new_pattern)
        if invalid_curies:
            entity2invalid_curies[entities.getId()] = invalid_curies

    else:
        for entity in tqdm(entities):
            invalid_curies = improve_uri_per_entity(entity, new_pattern)
            if invalid_curies:
                entity2invalid_curies[entity.getId()] = invalid_curies

            if type(entity) == UnitDefinition:
                for (
                    unit
                ) in (
                    entity.getListOfUnits()
                ):  # Unit needs to be handled within ListOfUnitDefinition
                    invalid_curies = improve_uri_per_entity(unit, new_pattern)
                    if invalid_curies:
                        entity2invalid_curies[unit.getId()] = invalid_curies

    return entity2invalid_curies


def polish_annotations(
    model: libModel, new_pattern: bool, outpath: str = None
) -> libModel:
    """Polishes all annotations in a model such that no duplicates are present
    & the same pattern is used for all CURIEs

    Args:
        - model (libModel):
            Model loaded with libSBML
        - new_pattern (bool):
            True if new pattern is wanted, otherwise False.
            Note that bioregistry internally only uses the new patter.
        - outpath (str, optional):
             Path to output file for invalid CURIEs detected by improve_uris
             Defaults to None.

    Returns:
        libModel:
            libSBML model with polished annotations
    """
    list_of_entity2invalid_curies = []
    listOf_dict = {
        "model": model,
        "compartment": model.getListOfCompartments(),
        "metabolite": model.getListOfSpecies(),
        "parameter": model.getListOfParameters(),
        "reaction": model.getListOfReactions(),
        "unit definition": model.getListOfUnitDefinitions(),
    }

    if model.isPackageEnabled("fbc"):
        listOf_dict["gene product"] = model.getPlugin("fbc").getListOfGeneProducts()

    if model.isPackageEnabled("groups"):
        listOf_dict["group"] = model.getPlugin("groups").getListOfGroups()

    # Adjust annotations in model
    for listOf in listOf_dict:
        print(f"Polish {listOf} annotations...")
        entity2invalid_curies = improve_uris(listOf_dict[listOf], new_pattern)
        list_of_entity2invalid_curies.append(entity2invalid_curies)

    all_entity2invalid_curies = reduce(
        lambda d1, d2: {**d1, **d2}, list_of_entity2invalid_curies
    )

    # Write invalid CURIEs to file if present
    if all_entity2invalid_curies:
        filename = (
            f'{model.getId()}_invalid_curies_{str(date.today().strftime("%Y%m%d"))}.csv'
        )
        if outpath:
            filename = Path(outpath, filename)
        else:
            filename = Path(filename)
        logging.warning(
            f"In the provided model {model.getId()} for {len(all_entity2invalid_curies)} entities invalid CURIEs were detected. "
            + f"These invalid CURIEs are saved to {filename}"
        )
        invalid_curies_df = parse_dict_to_dataframe(all_entity2invalid_curies)
        invalid_curies_df.columns = ["entity", "invalid_curie"]
        invalid_curies_df[["prefix", "identifier"]] = (
            invalid_curies_df.invalid_curie.str.split(r":", n=1, expand=True)
        )  # Required for identifiers that also contain a ':'
        invalid_curies_df = invalid_curies_df.drop("invalid_curie", axis=1)
        invalid_curies_df.to_csv(filename, index=False)

    return model


# Correct CVTerm qualifiers & qualifier types
# -------------------------------------------
def change_qualifier_per_entity(
    entity: SBase, new_qt, new_b_m_qt, specific_db_prefix: str = None
) -> list:
    """Updates Qualifiers to be MIRIAM compliant for an entity

    Args:
        - entity (SBase):
            A libSBML SBase object like model, GeneProduct, etc.
        - new_qt (Qualifier):
            A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - new_b_m_qt (QualifierType):
            A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - specific_db_prefix (str):
            Has to be set if only for a specific database the qualifier type should be changed. Can be 'kegg.genes', 'biocyc', etc.

    Returns:
        list:
            CURIEs that are not MIRIAM compliant
    """
    not_miriam_compliant = []
    pattern = f"{MIRIAM}|{OLD_MIRIAM}"
    cvterms = entity.getCVTerms()

    # for i in range(len(cvterms)):
    for cvterm in cvterms:
        tmp_set = SortedSet()
        sbo_set = SortedSet()

        if (
            cvterm.getBiologicalQualifierType() == 9
        ):  # 9 = BQB_OCCURS_IN (Reaction), Check for reactions with occursIn
            logging.info(
                f"CVTerm for {Fore.LIGHTYELLOW_EX}{str(entity)}{Style.RESET_ALL}"
                + f" is left as {Fore.LIGHTYELLOW_EX}{BiolQualifierType_toString(cvterm.getBiologicalQualifierType())}{Style.RESET_ALL}"
            )

        elif (
            cvterm.getModelQualifierType() == 1
        ):  # 1 = BQM_IS_DESCRIBED_BY (UnitDefinition), Check for UnitDefinitions with isDescribedBy
            logging.info(
                f"CVTerm for {Fore.LIGHTYELLOW_EX}{str(entity)}{Style.RESET_ALL}"
                + f" is left as {Fore.LIGHTYELLOW_EX}{ModelQualifierType_toString(cvterm.getModelQualifierType())}{Style.RESET_ALL}"
            )

        else:
            current_curies = [
                cvterm.getResourceURI(j) for j in range(cvterm.getNumResources())
            ]

            for cc in current_curies:

                current_curie = None

                if (specific_db_prefix != None) and (specific_db_prefix != ""):
                    if specific_db_prefix in cc:
                        current_curie = cc
                else:
                    current_curie = cc

                if (current_curie) and re.match(
                    pattern, current_curie, re.IGNORECASE
                ):  # If model contains identifiers without MIRIAM/OLD_MIRIAM these are kept
                    if re.search(r"sbo:", current_curie, re.IGNORECASE):
                        sbo_set.add(current_curie)
                    else:
                        tmp_set.add(current_curie)
                    cvterm.removeResource(current_curie)
                else:
                    not_miriam_compliant.append(current_curie)

            if sbo_set:
                add_uri_set(entity, BIOLOGICAL_QUALIFIER, BQB_HAS_PROPERTY, sbo_set)
            add_uri_set(entity, new_qt, new_b_m_qt, tmp_set)
            # cvterms.remove(i)

    if not_miriam_compliant:
        return not_miriam_compliant


def change_qualifiers(
    model: libModel,
    entity_type: str,
    new_qt,
    new_b_m_qt,
    specific_db_prefix: str = None,
) -> libModel:
    """Updates Qualifiers to be MIRIAM compliant for an entity type of a given model

    Args:
        - model (libModel):
            Model loaded with libSBML
        - entity_type (str):
            Any string of the following: model|compartment|metabolite|parameter|reaction|unit definition|unit|gene product|group
        - new_qt (Qualifier):
            A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - new_b_m_qt (QualifierType):
            A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - specific_db_prefix (str):
            Has to be set if only for a specific database the qualifier type should be changed. Can be 'kegg.genes', 'biocyc', etc.

    Returns:
        libModel:
            Model with changed qualifier for given entity type
    """
    not_miriam_compliant = []
    listOf_dict = {
        "model": model,
        "compartment": model.getListOfCompartments(),
        "metabolite": model.getListOfSpecies(),
        "parameter": model.getListOfParameters(),
        "reaction": model.getListOfReactions(),
        "unit definition": model.getListOfUnitDefinitions(),
    }

    if model.isPackageEnabled("fbc"):
        listOf_dict["gene product"] = model.getPlugin("fbc").getListOfGeneProducts()

    if model.isPackageEnabled("groups"):
        listOf_dict["group"] = model.getPlugin("groups").getListOfGroups()

    if entity_type == "model":  # Model needs to be handled like entity!
        not_miriam_compliant = change_qualifier_per_entity(
            listOf_dict.get("model"), new_qt, new_b_m_qt, specific_db_prefix
        )

    elif entity_type == "unit":
        for unit in listOf_dict.get(
            "unit definition"
        ):  # Unit needs to be handled within ListOfUnitDefinition
            not_miriam_compliant = change_qualifier_per_entity(
                unit, new_qt, new_b_m_qt, specific_db_prefix
            )

    else:
        try:
            for entity in tqdm(listOf_dict.get(entity_type)):
                not_miriam_compliant = change_qualifier_per_entity(
                    entity, new_qt, new_b_m_qt, specific_db_prefix
                )
        except TypeError:
            logging.info(
                "The entity " + entity_type + " is not present in " + model.getId()
            )

    if not_miriam_compliant:
        logging.warning(
            f"The following {len(not_miriam_compliant)} entities are not MIRIAM compliant: {not_miriam_compliant}"
        )

    return model


def change_all_qualifiers(model: libModel, lab_strain: bool) -> libModel:
    """Wrapper function to change qualifiers of all entities at once

    Args:
        - model (libModel):
            Model loaded with libSBML
        - lab_strain (bool):
            True if the strain was sequenced in a local lab

    Returns:
        libModel:
            Model with all qualifiers updated to be MIRIAM compliant
    """
    # Change all model entities to have the correct model qualifier
    entity_list_mod = ["model", "unit definition", "unit"]
    for entity in entity_list_mod:
        print(f"Change {str(entity)} qualifiers...")
        model = change_qualifiers(model, entity, MODEL_QUALIFIER, BQM_IS)

    # Change all remaining entities to have the correct biological qualifier
    entity_list = [
        "compartment",
        "metabolite",
        "parameter",
        "reaction",
    ]

    if model.isPackageEnabled("fbc"):
        entity_list.append("gene product")

    if model.isPackageEnabled("groups"):
        entity_list.append("group")

    for entity in entity_list:
        print(f"Change {str(entity)} qualifiers...")
        if lab_strain and entity == "gene product":
            model = change_qualifiers(
                model, "gene product", BIOLOGICAL_QUALIFIER, BQB_IS_HOMOLOG_TO
            )
        else:
            model = change_qualifiers(model, entity, BIOLOGICAL_QUALIFIER, BQB_IS)

    return model
