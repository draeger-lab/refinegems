#!/usr/bin/env python
"""General functions for curating a model

This module provides functionalities for curating models, including special functions for CarveMe models.

Since CarveMe version 1.5.1, the draft models from CarveMe contain pieces of information that are not correctly added to the
annotations. To address this, this module includes the following functionalities:

    - Add URIs from the entity IDs to the annotation field for metabolites & reactions
    - Transfer URIs from the notes field to the annotations for metabolites & reactions
    - Add URIs from the GeneProduct IDs to the annotations

The functionalities for CarveMe models, along with some of the following further functionalities, are gathered in the
main function :py:func:`~refinegems.curation.curate.polish_model`.

Further functionalities:

        - Setting boundary condition & constant for metabolites & reactions
        - Unit handling to add units & UnitDefinitions & to set units for parameters
        - Addition of default settings for compartments & metabolites
        - Addition of URIs to GeneProducts

            - via a mapping from model IDs to valid database IDs
            - via the KEGG API

        - Changing the CURIE pattern/CVTerm qualifier & qualifier type
        - Directionality control
"""

__author__ = "Famke Baeuerle and Carolin Brune and Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################

import cobra
import logging
import pandas as pd
import re
import warnings

from bioservices.kegg import KEGG
from cobra.io.sbml import _f_specie, _f_reaction
from libsbml import Model as libModel
from libsbml import GeneProduct, Species, ListOfSpecies, ListOfReactions, UnitDefinition
from pathlib import Path
from tqdm.auto import tqdm
from typing import Literal, Union

from .miriam import polish_annotations, change_all_qualifiers

from ..utility.cvterms import (
    add_cv_term_genes,
    add_cv_term_metabolites,
    add_cv_term_reactions,
    DB2PREFIX_METABS,
    DB2PREFIX_REACS,
    get_id_from_cv_term,
)
from ..utility.entities import get_gpid_mapping, create_fba_units, print_UnitDefinitions
from ..utility.io import load_a_table_from_database
from ..utility.util import DB2REGEX, test_biomass_presence

################################################################################
# variables
################################################################################

NH_PATTERN = re.compile(r"nh[3-4]")  #: :meta:

################################################################################
# functions
################################################################################

# sychronise annotations
# ----------------------


def update_annotations_from_others(model: libModel) -> libModel:
    """Synchronizes metabolite annotations for core, periplasm and extracelullar

    Args:
        - model (libModel):
            Model loaded with libSBML

    Returns:
        libModel:
            Modified model with synchronized annotations
    """
    for metab in model.getListOfSpecies():
        base = metab.getId()[:-2]
        for comp in ["_c", "_e", "_p"]:
            other_metab = model.getSpecies(base + comp)
            if other_metab is not None:
                if not other_metab.isSetMetaId():
                    other_metab.setMetaId("meta_" + other_metab.getId())
                for db_id, code in DB2PREFIX_METABS.items():
                    id = get_id_from_cv_term(metab, code)
                    for entry in id:
                        if entry is not None:
                            add_cv_term_metabolites(entry, db_id, other_metab)
    return model


def extend_gp_annots_via_mapping_table(
    model: libModel,
    mapping_tbl_file: Union[str,Path] = None,
    gff_paths: list[str] = None,
    email: str = None,
    contains_locus_tags: bool = False,
    lab_strain: bool = False,
    outpath: str = None,
) -> libModel:
    """
    | Extend GenePoduct annotations via mapping table.
    | If no mapping table is provided, a mapping table will be generated.

    Args:
        - model (libModel):
            Model loaded with libSBML
        - mapping_tbl_file (str|Path, optional):
            Path to a file containing a mapping table with columns ``model_id | X...`` where X can be ``REFSEQ``,
            ``NCBI``, ``locus_tag`` or ``UNCLASSIFIED``.
            The table can contain all of the ``X`` columns or at least one of them.
            Defaults to None.
        - gff_paths (list[str], optional):
            Path(s) to GFF file(s). Allowed GFF formats are: RefSeq, NCBI and Prokka.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - email (str, optional):
            E-mail for NCBI queries.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - contains_locus_tags (bool, optional):
            Specifies if provided model has locus tags within the label tag if set to True.
            This is only used when mapping_tbl_file == None.
            Defaults to False.
        - lab_strain (bool, optional):
            Specifies if a strain from no database was provided and thus has only homolog mappings if set to True.
            Defaults to False.
        - outpath (str, optional):
            Output path for location where the generated mapping table should be written to.
            This is only used when mapping_tbl_file == None.
            Defaults to None.

    Returns:
        libModel:
            Modified model with extended annotations for the GeneProducts
    """

    # 1. Get mapping
    # If no mapping table provided, get via function
    print("Get mapping information...")
    if not mapping_tbl_file:
        mapping_table = get_gpid_mapping(
            model, gff_paths, email, contains_locus_tags, outpath
        )
    else:  # Otherwise read in table from file
        mapping_tbl_file = mapping_tbl_file if isinstance(mapping_tbl_file, str) else str(mapping_tbl_file)
        mapping_table = pd.read_csv(mapping_tbl_file)

    # Drop all rows without model_id entries
    mapping_table.dropna(subset="model_id", inplace=True)
    # Use model_id as index
    mapping_table.set_index("model_id", inplace=True)
    # Replace all NaN values with empty string
    mapping_table.fillna("", inplace=True)

    # 2. Use table to fill in information in model
    # Get gene list
    gene_list = model.getPlugin("fbc").getListOfGeneProducts()

    print("Extending GeneProduct information...")
    for gene in tqdm(gene_list):
        # Get row of mapping table for current model_id
        gp_infos = mapping_table.loc[gene.getId(), :]

        # Add infos to current GeneProduct
        if gp_infos["name"]:
            gene.setName(gp_infos["name"])
        if gp_infos["locus_tag"]:
            gene.setLabel(gp_infos["locus_tag"])
        if ("REFSEQ" in gp_infos.index.to_list()) and gp_infos["REFSEQ"]:
            add_cv_term_genes(gp_infos["REFSEQ"], "REFSEQ", gene, lab_strain)
        if ("NCBI" in gp_infos.index.to_list()) and gp_infos["NCBI"]:
            add_cv_term_genes(gp_infos["NCBI"], "NCBI", gene, lab_strain)

    return model


def extend_gp_annots_via_KEGG(gene_list: list[GeneProduct], kegg_organism_id: str):
    """Adds KEGG gene & UniProt identifiers to the GeneProduct annotations

    Args:
        gene_list (list[GeneProduct]):
            libSBML ListOfGenes
        kegg_organism_id (str):
            Organism identifier in the KEGG database
    """
    k = KEGG()
    mapping_kegg_uniprot = k.conv("uniprot", kegg_organism_id)
    no_valid_kegg = []

    for gp in tqdm(gene_list):

        if gp.getId() != "G_spontaneous":
            kegg_gene_id = f"{kegg_organism_id}:{gp.getLabel()}"

            try:
                uniprot_id = mapping_kegg_uniprot[kegg_gene_id]

                add_cv_term_genes(kegg_gene_id, "KEGG", gp)
                add_cv_term_genes(uniprot_id.split(r"up:")[1], "UNIPROT", gp)

            except KeyError:
                no_valid_kegg.append(gp.getLabel())

    if no_valid_kegg:
        logging.info(
            f"The following {len(no_valid_kegg)} locus tags form no valid KEGG Gene ID: {no_valid_kegg} with the provided KEGG Organism ID: {kegg_organism_id}."
        )


def extend_metab_reac_annots_via_id(
    entity_list: Union[ListOfSpecies, ListOfReactions], id_db: str
) -> None:
    """Extends metabolite or reaction annotations by extracting the core ID from the entity ID and adding this ID as
    valid annotation if possible

    Args:
        - entity_list (Union[ListOfSpecies, ListOfReactions]):
            List of entities, either species (metabolites, ListOfSpecies) or reactions (ListOfReactions)
        - id_db (str):
            The database prefix to validate IDs against. Must correspond to a valid
            prefix in the Bioregistry

    Raises:
        - TypeError: Unsupported type for entity_list
    """
    # Initialise id_db correctly
    id_db = id_db.upper()

    # Set-up case-dependent variables
    match entity_list:
        case ListOfSpecies():
            bigg_db = "bigg_metabolites"
            bigg_id_type = "universal_bigg_id"
            add_cv_term = add_cv_term_metabolites
            id_handler = _f_specie
            db2prefix = DB2PREFIX_METABS
        case ListOfReactions():
            bigg_db = "bigg_reactions"
            bigg_id_type = "id"
            add_cv_term = add_cv_term_reactions
            id_handler = _f_reaction
            db2prefix = DB2PREFIX_REACS
        case _:
            raise TypeError(
                f"Unsupported type for entity_list {type(entity_list)}. Must be ListOfSpecies or ListOfReactions."
            )

    # Set-up default variables
    try:
        id_db_prefix = db2prefix[id_db]
    except KeyError:
        print(
            f"""
KeyError: with id_db=\'{id_db}\'
id_db must be one of the valid database names: {db2prefix.keys()}
If your id_db is not part of the list, please contact the developers.
              """
        )
        return
    try:
        db_pattern = DB2REGEX[id_db_prefix]
    except KeyError:
        print(
            f"KeyError: id_db_prefix = {id_db_prefix} must be one of the valid prefixes in https://bioregistry.io/."
        )
        return

    # Get BiGG IDs for VMH ID == BiGG ID validation
    if "vmh" in id_db_prefix.lower():
        bigg_ids = load_a_table_from_database(f"SELECT bigg_id FROM {bigg_db}")
        bigg_ids = set(bigg_ids[bigg_id_type].tolist())

    # Get ID from entity & add annotation if valid database ID
    for entity in entity_list:

        # Get ID
        current_id = entity.getId()

        # Use current_id as metaid if no metaid is present
        if not entity.isSetMetaId():
            entity.setMetaId(f"meta_{current_id}")

        if re.fullmatch(db_pattern, current_id, re.IGNORECASE):

            # Remove prefix
            current_id = id_handler(current_id)

            # Unset annotations if no CV terms exist
            if entity.getNumCVTerms() == 0:
                entity.unsetAnnotation()

            # Get ID for annotation
            if isinstance(entity, Species):  # Remove compartment suffix
                id_for_anno = re.sub(
                    f"(_|\[){entity.getCompartment()}\]?$", "", current_id
                )
            else:
                id_for_anno = current_id

            # Add ID as URI to annotation
            add_cv_term(id_for_anno, id_db, entity)

            # Add BiGG ID to annotation, additionally, if VMH and valid BiGG ID
            if ("vmh" in id_db_prefix.lower()) and (id_for_anno in bigg_ids):
                add_cv_term(id_for_anno, id_db, entity)


def extend_metab_reac_annots_via_notes(
    entity_list: Union[ListOfSpecies, ListOfReactions],
) -> None:
    """Extends metabolite or reaction annotations by extracting valid URIs from the notes section of the provided
    entities and cleans up the notes by removing processed elements

    Args:
        - entity_list (Union[ListOfSpecies, ListOfReactions]):
            List of entities, either species (metabolites, ListOfSpecies) or reactions (ListOfReactions)
    Raises:
        - TypeError: Unsupported type for entity_list
    """

    # Set-up case-dependent variables
    match entity_list:
        case ListOfSpecies():
            db2prefix = DB2PREFIX_METABS
            add_cv_term = add_cv_term_metabolites
        case ListOfReactions():
            db2prefix = DB2PREFIX_REACS
            add_cv_term = add_cv_term_reactions
        case _:
            raise TypeError(
                f"Unsupported type for entity_list {type(entity_list)}. Must be ListOfSpecies or ListOfReactions."
            )

    # Get notes & add annotation from notes to CVTerms if valid URI
    for entity in entity_list:

        if not entity.isSetMetaId():
            entity.setMetaId("meta_" + entity.getId())

        notes_list = []
        elem_used = []
        notes_string = entity.getNotesString().split("\n")
        for elem in notes_string:
            for id_db in db2prefix.keys():
                if "<p>" + id_db in elem.upper():
                    elem_used.append(elem)
                    # @DEBUG print(elem.strip()[:-4].split(r': ')[1])
                    fill_in = re.split(r":\s*", elem.strip()[:-4])[1]
                    if (";") in fill_in:
                        if not re.search(r"inchi", id_db, re.IGNORECASE):
                            entries = fill_in.split(r";")
                            for entry in entries:
                                if not re.fullmatch(
                                    r"^nan$", entry.strip(), re.IGNORECASE
                                ):
                                    add_cv_term(entry.strip(), id_db, entity)
                    else:
                        if not re.fullmatch(r"^nan$", fill_in, re.IGNORECASE):
                            add_cv_term(fill_in.strip(), id_db, entity)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = " ".join([str(elem) + "\n" for elem in notes_list])
        entity.unsetNotes()
        entity.setNotes(new_notes)
        # @DEBUG print(species.getAnnotationString())


# correct basic model set-up
# ---------------------------


def polish_model_units(model: libModel) -> None:
    """Replaces the list of unit definitions with the unit definitions needed for FBA:

       - mmol per gDW per h
       - mmol per gDW
       - hour (h)
       - femto litre (fL)

    Args:
        - model (libModel):
            Model loaded with libSBML
    """

    # Get FBA unit definitions per refineGEMs definition
    fba_unit_defs = create_fba_units(model)

    # Get model unit definitions
    model_unit_defs = model.getListOfUnitDefinitions().clone()

    # If list of unit definitions is not empty, replace all units with the defined FBA units
    # & Print the non-FBA unit definitions
    if model_unit_defs:

        # List to collect all non-FBA unit definitions
        removed_unit_defs = []

        # Check if model unit definitions fit to the fba unit definitions
        for model_ud in model_unit_defs:
            for fba_ud in fba_unit_defs:

                # In case of identical unit definitions, remove unit definition from fba unit def list
                if not UnitDefinition.areIdentical(fba_ud, model_ud):
                    removed_unit_defs.append(model_ud)

        # Remove all model unit definitions
        model.getListOfUnitDefinitions().clear(doDelete=True)

        # Only print list if UnitDefinitions were removed
        if len(removed_unit_defs) == 0:
            logging.warning(
                """
            The following UnitDefinition objects were removed. 
            The reasoning is that
            \t(a) these UnitDefinitions are not contained in the UnitDefinition list of this program and
            \t(b) the UnitDefinitions defined within this program are handled as ground truth.
            Thus, the following UnitDefinitions are not seen as relevant for the model.
            """
            )
            print_UnitDefinitions(removed_unit_defs)

    # Add all defined FBA units to the model
    for unit_def in fba_unit_defs:
        model.getListOfUnitDefinitions().append(unit_def)


def set_model_default_units(model: libModel):
    """Sets default units of model

    Args:
        - model (libModel):
            Model loaded with libSBML
    """
    for unit in model.getListOfUnitDefinitions():

        unit_id = unit.getId()

        if re.fullmatch(r"mmol_per_gDW", unit_id, re.IGNORECASE):

            if not (model.isSetExtentUnits() and model.getExtentUnits() == unit_id):
                model.setExtentUnits(unit_id)

            if not (
                model.isSetSubstanceUnits() and model.getSubstanceUnits() == unit_id
            ):
                model.setSubstanceUnits(unit_id)

        if not (
            model.isSetTimeUnits() and model.getTimeUnits() == unit_id
        ) and re.fullmatch(r"hr?", unit_id, re.IGNORECASE):
            model.setTimeUnits(unit_id)

        if not (
            model.isSetVolumeUnits() and model.getVolumeUnits() == unit_id
        ) and re.fullmatch(r"fL", unit_id, re.IGNORECASE):
            model.setVolumeUnits(unit_id)


def set_units_of_parameters(model: libModel):
    """Sets units of parameters in model

    Args:
        - model (libModel):
            Model loaded with libSBML
    """
    for (
        param
    ) in (
        model.getListOfParameters()
    ):  # needs to be added to list of unit definitions aswell
        if any(
            (
                unit_id := re.fullmatch(
                    r"mmol_per_gDW_per_hr?", unit.getId(), re.IGNORECASE
                )
            )
            for unit in model.getListOfUnitDefinitions()
        ):
            if not (param.isSetUnits() and param.getUnits() == unit_id.group(0)):
                param.setUnits(unit_id.group(0))


def add_compartment_structure_specs(model: libModel):
    """| Adds the required specifications for the compartment structure
    | if not set (size & spatial dimension)

    Args:
        - model (libModel):
            Model loaded with libSBML
    """
    for compartment in model.getListOfCompartments():

        if not compartment.isSetSize():
            compartment.setSize(float("NaN"))

        if not compartment.isSetSpatialDimensions():
            compartment.setSpatialDimensions(3)

        if any(
            (unit_id := re.fullmatch(r"fL", unit.getId(), re.IGNORECASE))
            for unit in model.getListOfUnitDefinitions()
        ):
            if not (
                compartment.isSetUnits() and compartment.getUnits() == unit_id.group(0)
            ):
                compartment.setUnits(unit_id.group(0))


def set_initial_amount_metabs(model: libModel):
    """Sets initial amount to all metabolites if not already set or if initial concentration is not set

    Args:
        - model (libModel):
            Model loaded with libSBML
    """
    for species in model.getListOfSpecies():

        if not (species.isSetInitialAmount() or species.isSetInitialConcentration()):
            species.setInitialAmount(float("NaN"))


def polish_entity_conditions(entity_list: Union[ListOfSpecies, ListOfReactions]):
    """Sets boundary condition and constant if not set for an entity

    Args:
        - entity_list (Union[ListOfSpecies, ListOfReactions]):
            libSBML ListOfSpecies or ListOfReactions
    """
    match entity_list:
        case ListOfSpecies():
            for entity in entity_list:
                if not entity.getBoundaryCondition():
                    entity.setBoundaryCondition(False)
                if not entity.getConstant():
                    entity.setConstant(False)
        case ListOfReactions():
            pass
        case _:
            logging.warning(
                f"Unsupported type for entity_list {type(entity_list)}. Must be ListOfSpecies or ListOfReactions."
            )


# duplicates
# ----------


def resolve_duplicate_reactions(
    model: cobra.Model, based_on: str = "reaction", remove_reac: bool = True
) -> cobra.Model:
    """Resolve and remove duplicate reaction based on their reaction equation
    and matching database identifiers. Only if all match or a comparison with nan occurs will one of
    the reactions be removed.

    Args:
        - model (cobra.Model):
            A model loaded with COBRApy.
        - based_on (str, optional):
            Label to base the resolvement process on .
            Can be 'reaction' or any other annotation label.
            Defaults to 'reaction'.
        - remove_reac (bool, optional):
            When True, combines and remove duplicates.
            Otherwise only reports the findings.
            Defaults to True.

    Returns:
        cobra.Model:
            The model.
    """

    # get annotation and compartment information
    anno_reac = []
    for r in model.reactions:
        anno_reac.append(
            {"id": r.id, "compartment": str(r.compartments), "reaction": r.reaction}
            | r.annotation
        )
    df_reac = pd.DataFrame.from_dict(anno_reac)

    # check if based_on is valid
    if not based_on in df_reac.columns.tolist():
        warnings.warn(
            f"Warning: Annotation column {based_on} does not exists. Search for duplicates will be skipped."
        )
        return model

    # set basic parameters
    skip_cols = ["id", "compartment", "bigg.reaction", "reaction", based_on]
    colnames = df_reac.columns.tolist()

    for c in df_reac.groupby("compartment"):
        # note: using groupby drops nans
        for mnx in c[1].groupby(based_on):
            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                # check annotations
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = (
                            mnx[1][col].dropna().explode().unique().tolist()
                        )
                    else:
                        dupl = False
                        break

                # if duplicate found
                if dupl:
                    if remove_reac:

                        # choose reaction to keep
                        keep_reac = model.reactions.get_by_id(mnx[1]["id"].tolist()[0])

                        # resolve annotations
                        for key, value in annotations.items():
                            if len(value) > 0 and not key in keep_reac.annotation:
                                keep_reac.annotation[key] = value

                        # combine gene reaction rules
                        for r_id in mnx[1]["id"].tolist()[1:]:
                            gpr_to_add = model.reactions.get_by_id(
                                r_id
                            ).gene_reaction_rule
                            if gpr_to_add and gpr_to_add != "":
                                if (
                                    keep_reac.gene_reaction_rule
                                    and keep_reac.gene_reaction_rule != ""
                                ):  # add two existing rules
                                    keep_reac.gene_reaction_rule = (
                                        keep_reac.gene_reaction_rule
                                        + " or "
                                        + gpr_to_add
                                    )
                                else:
                                    keep_reac.gene_reaction_rule = (
                                        gpr_to_add  # add the one existing the other
                                    )
                            else:
                                pass  # nothing to add
                            model.reactions.get_by_id(r_id).remove_from_model()
                            print(
                                f"\tDuplicate reaction {r_id} found. Combined to {keep_reac.id} and deleted."
                            )
                    else:
                        print(
                            f'\tDuplicate reactions {", ".join(mnx[1]["id"].tolist())} found.'
                        )

    return model


def resolve_duplicate_metabolites(
    model: cobra.Model, based_on: str = "metanetx.chemical", replace: bool = True
) -> cobra.Model:
    """Resolve duplicate metabolites in a model. Metabolites are considered
    duplicate if they share the same annotations (same or nan).

    .. note::
        Depending on the starting database, the results might differ.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - based_on (str, optional):
            Label to base the resolvement process on .
            Can be any annotation label.
            Defaults to 'metanetx.chemical'.
        - replace (bool, optional):
            Either report the duplicates (False)
            or replace them with one (True). Defaults to True.

    Returns:
        cobra.Model:
            The model.
    """

    # get annotation and compartment information
    anno_meta = []
    for m in model.metabolites:
        anno_meta.append({"id": m.id, "compartment": m.compartment} | m.annotation)
    df_meta = pd.DataFrame.from_dict(anno_meta)

    # check if based_on is valid
    if not based_on in df_meta.columns.tolist():
        warnings.warn(
            f"Warning: Annotation {based_on} not found. Search for metabolite duplicates skipped."
        )
        return model

    # set basic parameters
    skip_cols = ["id", "compartment", "bigg.metabolite", based_on]
    colnames = df_meta.columns.tolist()
    # get objective function
    bof_list = test_biomass_presence(model)
    if len(bof_list) == 1:
        objective_function = bof_list[0]
    elif len(bof_list) > 1:
        mes = f"Multiple BOFs detected. Will be using {bof_list[0]}"
        warnings.warn(mes, category=UserWarning)
        objective_function = bof_list[0]
    else:
        mes = "No BOF detected. Might lead to problems during duplicate removal."
        warnings.warn(mes, category=UserWarning)

    for c in df_meta.groupby("compartment"):
        # note: using groupby drops nans
        # note: mnx as starting point was choosen, as it seems to have the most annotations (easy to get)
        c[1][based_on] = c[1][based_on].apply(
            lambda x: tuple(x) if isinstance(x, list) else x
        )
        for mnx in c[1].groupby(based_on):

            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    # check annotation, if current group truly consists of duplicates
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = (
                            mnx[1][col].dropna().explode().unique().tolist()
                        )
                    else:
                        dupl = False
                        break
            # no duplicates = check next
            else:
                continue

            # if duplicate is found:
            if dupl:
                # either replace ...
                if replace:
                    # choose metabolite to keep
                    keep_meta = model.metabolites.get_by_id(mnx[1]["id"].tolist()[0])
                    # resolve annotations
                    for key, value in annotations.items():
                        if len(value) > 0 and not key in keep_meta.annotation:
                            keep_meta.annotations[key] = value
                    # note: charge and formula should be in valid range and be corrected by MCC (if needed)

                    # replace duplicates with the metabolite to be kept
                    #     to ensure consistency, only delete duplicate metabolites, which
                    #     do NOT share ANY reactions

                    # retrieve reactions for metabolites set for keeping
                    keep_reac = [
                        _.id
                        for _ in model.metabolites.get_by_id(keep_meta.id).reactions
                    ]
                    # iterate over metabolites set for deletion
                    for del_meta_id in mnx[1]["id"].tolist()[1:]:
                        # retrieve reaction for metabolites set for deletion
                        del_reac = [
                            _.id
                            for _ in model.metabolites.get_by_id(del_meta_id).reactions
                        ]
                        # get intersection of reactions (keep + del)
                        reac_intersec = list(set(keep_reac) & set(del_reac))

                        # if intersection empty, metabolite is with a high probability indeed a duplicate
                        # Special case: NH3 / NH4
                        #    intersection does not have to be emtpy, know 'problem' caused by CarveMe
                        if len(reac_intersec) == 0 or all(
                            [
                                re.search(NH_PATTERN, _)
                                for _ in [keep_meta.id, del_meta_id]
                            ]
                        ):
                            # automated deletion is only advisable, if consistency can be retained
                            perform_deletion = True
                            with model as model_del:

                                # if the special case is detected ...
                                if all(
                                    [
                                        re.search(NH_PATTERN, _)
                                        for _ in [keep_meta.id, del_meta_id]
                                    ]
                                ):
                                    print(
                                        f"\tSpecial case -Duplicate NH4/NH3- detected.\n\tTrying to solve by additionally removing reactions containing both metabolites."
                                    )
                                    # ... remove reactions with nh3 and nh4 both present
                                    for del_reac_id in reac_intersec:
                                        # if objective_function is part of the set
                                        # automated deletion is (currently) not possible
                                        if del_reac_id == objective_function:
                                            perform_deletion = False
                                            break
                                        model_del.reactions.get_by_id(
                                            del_reac_id
                                        ).remove_from_model()
                                    # set the metabolites to be deleted to be the one NOT in the objective functions
                                    # to avoid inconsistencies
                                    if del_meta_id in [
                                        _.id
                                        for _ in model_del.reactions.get_by_id(
                                            objective_function
                                        ).metabolites
                                    ]:
                                        temp = del_meta_id
                                        del_meta_id = keep_meta.id
                                        keep_meta = model_del.metabolites.get_by_id(
                                            temp
                                        )

                                # try replacing metabolite with the kept duplicate ...
                                for reac in model_del.metabolites.get_by_id(
                                    del_meta_id
                                ).reactions:

                                    reac.reaction = reac.reaction.replace(
                                        del_meta_id, keep_meta.id
                                    )

                                    # skip if objective function is found
                                    if reac.id == objective_function:
                                        continue

                                    # check if consistency is still intact
                                    balance_test = reac.check_mass_balance()
                                    if not reac.boundary and len(balance_test) > 0:
                                        # try fixing H-balance
                                        if "H" in balance_test.keys():
                                            # ..............................
                                            # TODO:
                                            #    get H according to compartment
                                            #    current implementation relies heavily
                                            #    on 'correct' use input: compartment should have format C_c or c (C_p, p, C_e, e etc.)
                                            # ..............................
                                            reac_comp = reac.compartments.pop()[-1]
                                            if reac_comp == "c":
                                                reac.subtract_metabolites(
                                                    {"h_c": balance_test["H"]}
                                                )
                                            elif reac_comp == "p":
                                                reac.subtract_metabolites(
                                                    {"h_p": balance_test["H"]}
                                                )
                                            elif reac_comp == "e":
                                                reac.subtract_metabolites(
                                                    {"h_e": balance_test["H"]}
                                                )
                                            else:
                                                perform_deletion = False
                                                break
                                        # ..............................
                                        # TODO:
                                        #    fix other possible problems
                                        # ..............................

                                        # finally, check balance again (continue only if fixed, else break)
                                        if len(reac.check_mass_balance()) > 0:
                                            perform_deletion = False
                                            break

                                    else:
                                        continue

                                # if not problems are found, duplicate is removed
                                if perform_deletion:
                                    model = model_del.copy()
                                    print(
                                        f"\tDuplicate metabolite {del_meta_id} found. Replaced with {keep_meta.id}."
                                    )
                                # if problems are not solvable, duplicate is kept and only reported
                                else:
                                    print(
                                        f"\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tAutomated deletion not possible due to problems with consistency."
                                    )

                        # else, metabolite is kept
                        #       since it might be an isomer, elongation, or other explanation
                        #       for the same annotation
                        else:
                            print(
                                f"\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tKept, as reaction containing both metabolites was found."
                            )

                # ... or only report duplicates
                else:
                    print(
                        f'\tDuplicate metabolite(s) {", ".join(mnx[1]["id"].tolist())} found.'
                    )

    return model


def resolve_duplicates(
    model: cobra.Model,
    check_reac: bool = True,
    check_meta: Literal["default", "exhaustive", "skip"] = "default",
    replace_dupl_meta: bool = True,
    remove_unused_meta: bool = False,
    remove_dupl_reac: bool = True,
) -> cobra.Model:
    """Resolve and remove (optional) duplicate metabolites and reactions in the model.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - check_reac (bool, optional):
            Whether to check reactions for duplicates.
            Defaults to True.
        - check_meta (Literal['default','exhaustive','skip'], optional):
            Whether to check for duplicate metabolites.
            Defaults to 'default'.
        - replace_dupl_meta (bool, optional):
            Option to replace/remove duplicate metabolites.
            Defaults to True.
        - remove_unused_meta (bool, optional):
            Option to remove unused metabolites.
            Defaults to False.
        - remove_dupl_reac (bool, optional):
            Option to combine/remove duplicate reactions.
            Defaults to True.

    Returns:
        cobra.Model:
            The (edited) model.
    """

    # resolve duplicate metabolites
    if check_meta == "default":
        # resolve duplicates starting with the metanetx.chemical database identifiers
        model = resolve_duplicate_metabolites(model, replace=replace_dupl_meta)
    elif check_meta == "exhaustive":
        # resolve duplicates by starting at every database identifer one after another
        # note: bigg and sbo are skipped as sbo gives not much information and bigg is
        #       usually the one that differs (naming issue)
        anno_types = set()
        # get all database annotation types present in the model
        for m in model.metabolites:
            anno_types = anno_types | set(m.annotation.keys())
        for colname in [_ for _ in anno_types if not _ in ["bigg.metabolite", "sbo"]]:
            model = resolve_duplicate_metabolites(
                model, colname, replace=replace_dupl_meta
            )
    elif check_meta == "skip":
        print("\tSkip check for duplicate metabolites.")
    else:
        warnings.warn(
            f"Warning: Unknown option for metabolites duplicate checking {check_meta}. Search for metabolite duplicates skipped."
        )

    # remove now unused metabolites
    if remove_unused_meta:
        model, removed = cobra.manipulation.delete.prune_unused_metabolites(model)
        print(
            f'\tThe following metabolites () have been removed: {", ".join([x.id for x in removed])}'
        )

    # resolve duplicate reactions
    if check_reac:
        model = resolve_duplicate_reactions(
            model, based_on="reaction", remove_reac=remove_dupl_reac
        )

    return model


# Directionality Control
# ----------------------


def check_direction(model: cobra.Model, data: Union[pd.DataFrame, str]) -> cobra.Model:
    """Check the direction of reactions by searching for matching MetaCyc,
    KEGG and MetaNetX IDs as well as EC number in a downloaded BioCyc (MetaCyc)
    database table or dataFrame (need to contain at least the following columns:
    Reactions (MetaCyc ID),EC-Number,KEGG reaction,METANETX,Reaction-Direction.

    Args:
        model (cobra.Model):
            The model loaded with COBRApy.
        data (pd.DataFrame | str):
            Either a pandas DataFrame or a path to a CSV file
            containing the BioCyc smart table.

    Raises:
        - TypeError: Unknown data type for parameter data

    Returns:
        cobra.Model:
            The edited model.
    """

    match data:
        # already a DataFrame
        case pd.DataFrame():
            pass
        case str():
            # load from a table
            data = pd.read_csv(data, sep="\t")
            # rewrite the columns into a better comparable/searchable format
            data["KEGG reaction"] = data["KEGG reaction"].str.extract(r".*>(R\d*)<.*")
            data["METANETX"] = data["METANETX"].str.extract(r".*>(MNXR\d*)<.*")
            data["EC-Number"] = data["EC-Number"].str.extract(r"EC-(.*)")
        case _:
            mes = f"Unknown data type for parameter data: {type(data)}"
            raise TypeError(mes)

    # check direction
    # --------------------
    for r in model.reactions:

        direction = None
        # easy case: metacyc is already (corretly) annotated
        if (
            "metacyc.reaction" in r.annotation
            and len(data[data["Reactions"] == r.annotation["metacyc.reaction"]]) != 0
        ):
            direction = data[data["Reactions"] == r.annotation["metacyc.reaction"]][
                "Reaction-Direction"
            ].iloc[0]
            r.notes["BioCyc direction check"] = f"found {direction}"
        # complicated case: no metacyc annotation
        else:
            annotations = []

            # collect matches
            if (
                "kegg.reaction" in r.annotation
                and r.annotation["kegg.reaction"] in data["KEGG reaction"].tolist()
            ):
                annotations.append(
                    data[data["KEGG reaction"] == r.annotation["kegg.reaction"]][
                        "Reactions"
                    ].tolist()
                )
            if (
                "metanetx.reaction" in r.annotation
                and r.annotation["metanetx.reaction"] in data["METANETX"].tolist()
            ):
                annotations.append(
                    data[data["METANETX"] == r.annotation["metanetx.reaction"]][
                        "Reactions"
                    ].tolist()
                )
            if (
                "ec-code" in r.annotation
                and r.annotation["ec-code"] in data["EC-Number"].tolist()
            ):
                annotations.append(
                    data[data["EC-Number"] == r.annotation["ec-code"]][
                        "Reactions"
                    ].tolist()
                )

            # check results
            # no matches
            if len(annotations) == 0:
                r.notes["BioCyc direction check"] = "not found"

            # matches found
            else:
                # built intersection
                intersec = set(annotations[0]).intersection(*annotations)
                # case 1: exactly one match remains
                if len(intersec) == 1:
                    entry = intersec.pop()
                    direction = data[data["Reactions"] == entry][
                        "Reaction-Direction"
                    ].iloc[0]
                    r.annotation["metacyc.reaction"] = entry
                    r.notes["BioCyc direction check"] = f"found {direction}"

                # case 2: multiple matches found -> inconclusive
                else:
                    r.notes["BioCyc direction check"] = f"found, but inconclusive"

        # update direction if possible and needed
        if not pd.isnull(direction):
            if "REVERSIBLE" in direction:
                # set reaction as reversible by setting default values for upper and lower bounds
                r.lower_bound = cobra.Configuration().lower_bound
            elif "RIGHT-TO-LEFT" in direction:
                # invert the default values for the boundaries
                r.lower_bound = cobra.Configuration().lower_bound
                r.upper_bound = 0.0
            elif "LEFT-To-RIGHT" in direction:
                # In case direction was already wrong
                r.lower_bound = 0.0
                r.upper_bound = cobra.Configuration().upper_bound
            else:
                # left to right case is the standard for adding reactions
                # = nothing left to do
                continue

    return model


# Perform all clean-up steps
# --------------------------
def polish_model(
    model: libModel,
    id_db: str = "BiGG",
    mapping_tbl_file: str = None,
    gff_paths: list[str] = None,
    email: str = None,
    contains_locus_tags: bool = False,
    lab_strain: bool = False,
    kegg_organism_id: str = None,
    reaction_direction: str = None,
    outpath: str = None,
) -> libModel:
    """Completes all steps to polish a model

    .. note::

        So far only tested for models having either BiGG or VMH identifiers.

    Args:
        - model (libModel):
             Model loaded with libSBML
        - id_db (str, optional):
            Main database where identifiers in model come from.
            Defaults to 'BiGG'.
        - mapping_tbl_file (str, optional):
            Path to a file containing a mapping table with columns ``model_id | X...`` where X can be ``REFSEQ``,
            ``NCBI``, ``locus_tag`` or ``UNCLASSIFIED``.
            The table can contain all of the ``X`` columns or at least one of them.
            Defaults to None.
        - gff_paths (list[str], optional):
            Path(s) to GFF file(s). Allowed GFF formats are: RefSeq, NCBI and Prokka.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - email (str, optional):
            E-mail for NCBI queries.
            This is only used when mapping_tbl_file == None.
            Defaults to None.
        - contains_locus_tags (bool, optional):
            Specifies if provided model has locus tags within the label tag if set to True.
            This is only used when mapping_tbl_file == None.
            Defaults to False.
        - lab_strain (bool, optional):
            Specifies if a strain from no database was provided and thus has only homolog mappings, if set to True.
            Defaults to False.
        - kegg_organism_id (str, optional):
            KEGG organism identifier if available.
            Defaults to None.
        - reaction_direction (str, optional):
            Path to a CSV file containing the BioCyc smart table with the columns
            ``Reactions (MetaCyc ID) | EC-Number | KEGG reaction | METANETX | Reaction-Direction``.
            For more details see :py:func:`~refinegems.curation.curate.check_direction`
            Defaults to None.
        - outpath (str, optional):
            Output path for mapping table from model ID to valid database IDs (if mapping_tbl_file == None)
            & incorrect annotations file(s).
            Defaults to None.

    Returns:
        libModel:
            Polished libSBML model
    """
    ### Set-up
    # Get ListOf objects
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    gene_list = model.getPlugin("fbc").getListOfGeneProducts()

    ### unit definition ###
    polish_model_units(model)
    set_model_default_units(model)
    set_units_of_parameters(model)
    add_compartment_structure_specs(model)
    set_initial_amount_metabs(model)

    ### improve metabolite, reaction and gene annotations ###
    extend_metab_reac_annots_via_id(metab_list, id_db)
    extend_metab_reac_annots_via_id(reac_list, id_db)
    extend_metab_reac_annots_via_notes(metab_list)
    extend_metab_reac_annots_via_notes(reac_list)
    update_annotations_from_others(model)

    ### Extend annotations for GeneProducts ###
    extend_gp_annots_via_mapping_table(
        model,
        mapping_tbl_file,
        gff_paths,
        email,
        contains_locus_tags,
        lab_strain,
        outpath
    )
    if kegg_organism_id:
        extend_gp_annots_via_KEGG(gene_list, kegg_organism_id)

    ### Check reaction direction ###
    if reaction_direction:
        check_direction(model, reaction_direction)

    ### set boundaries and constants ###
    polish_entity_conditions(metab_list)
    polish_entity_conditions(reac_list)

    ### MIRIAM compliance of CVTerms ###
    print(
        "Remove duplicates & transform all CURIEs to the new identifiers.org pattern (: between db and ID):"
    )
    model = polish_annotations(model, True, outpath)
    print("Changing all qualifiers to be MIRIAM compliant:")
    model = change_all_qualifiers(model, lab_strain)

    return model
