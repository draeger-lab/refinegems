#!/usr/bin/env python
"""Provides functions to investigate the model and test with MEMOTE

These functions enable simple testing of any model using MEMOTE and access to its number of reactions, metabolites and genes.
"""

__author__ = "Famke Baeuerle and Alina Renz and Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra

from libsbml import Model as libModel
from cobra import Model as cobraModel

from memote.support import consistency

# needed by memote.support.consistency
from memote.support import consistency_helpers as con_helpers

from ..utility.util import test_biomass_presence

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# get basic model info
# --------------------


def get_reac_with_gpr(model: cobra.Model) -> tuple[list[str], list[str]]:
    """Extract the reactions that have a gene production rule from a given model.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.

    Returns:
        tuple:
            Lists of reactions with GPR (1) & (2):

            (1) list: List of normal reactions with gpr
            (2) list: List of pseudoreactions with gpr
    """
    pseudoreactions = model.boundary
    for biomass in test_biomass_presence(model):
        pseudoreactions.append(model.reactions.get_by_id(biomass))

    normal = []
    pseudo = []
    for reac in model.reactions:
        # check for GPR
        if len(reac.genes) > 0:
            if reac in pseudoreactions:
                pseudo.append(reac.id)
            else:
                normal.append(reac.id)

    return (normal, pseudo)


def get_orphans_deadends_disconnected(
    model: cobraModel,
) -> tuple[list[str], list[str], list[str]]:
    """Uses MEMOTE functions to extract orphans, deadends and disconnected metabolites

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        tuple:
            Lists of metabolites that might cause errors (1) - (3)

            (1) list: List of orphans
            (2) list: List of deadends
            (3) list: List of disconnected metabolites
    """
    orphans = consistency.find_orphans(model)
    deadends = consistency.find_deadends(model)
    disconnected = consistency.find_disconnected(model)

    orphan_list = []
    if len(orphans) > 0:
        for orphan in orphans:
            orphan_list.append(orphan.id)

    deadend_list = []
    if len(deadends) > 0:
        for deadend in deadends:
            deadend_list.append(deadend.id)

    disconnected_list = []
    if len(disconnected) > 0:
        for disc in disconnected:
            disconnected_list.append(disc.id)

    return orphan_list, deadend_list, disconnected_list


def get_mass_charge_unbalanced(model: cobraModel) -> tuple[list[str], list[str]]:
    """Creates lists of mass and charge unbalanced reactions, without exchange reactions since they are unbalanced per definition

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        tuple:
            Lists of reactions that might cause errors (1) & (2)

            (1) list: List of mass unbalanced reactions
            (2) list: List of charge unbalanced reactions
    """

    mass_unbalanced = consistency.find_mass_unbalanced_reactions(model.reactions)
    charge_unbalanced = consistency.find_charge_unbalanced_reactions(model.reactions)

    mass_list = []
    if len(mass_unbalanced) > 0:
        for reac in mass_unbalanced:
            if reac.id not in model.exchanges:
                mass_list.append(reac.id)

    charge_list = []
    if len(charge_unbalanced) > 0:
        for reac in charge_unbalanced:
            if reac.id not in model.exchanges:
                charge_list.append(reac.id)

    return mass_list, charge_list


# other
# -----


def get_metabs_with_one_cvterm(model: libModel) -> list[str]:
    """Reports metabolites which have only one annotation, can be used as basis for further annotation research

    Args:
        - model (libModel):
            Model loaded with libSBML

    Returns:
        list:
            Metabolite Ids with only one annotation
    """
    spe = model.getListOfSpecies()

    only_one = []  # safe metab with only BiGG annotation
    for sb in spe:
        if sb.isSetId():
            pid = sb.getId()
            if sb.getCVTerm(0).getNumResources() == 1:
                only_one.append(pid)

    return only_one


# SBO terms
# ---------


def get_reactions_per_sbo(model: libModel) -> dict:
    """Counts number of reactions of all SBO Terms present

    Args:
        - model (libModel):
            Model loaded with libSBML

    Returns:
        dict:
            SBO Term as keys and number of reactions as values
    """
    sbos_dict = {}
    for react in model.getListOfReactions():
        sbo = react.getSBOTerm()
        if sbo in sbos_dict.keys():
            sbos_dict[sbo] += 1
        else:
            sbos_dict[sbo] = 1
    return sbos_dict
