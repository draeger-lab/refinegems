#!/usr/bin/env python
""" Provides functions to investigate the model and test with MEMOTE

These functions enable simple testing of any model using MEMOTE and access to its number of reactions, metabolites and genes.
"""

__author__ = "Famke Baeuerle and Alina Renz and Carolin Brune"

################################################################################
# requirements
################################################################################

import memote
import json
import pandas as pd
import numpy as np
import cobra
import time

from cobra import Reaction
from libsbml import Model as libModel
from cobra import Model as cobraModel
from typing import Literal

from memote.support import consistency
# needed by memote.support.consistency
from memote.support import consistency_helpers as con_helpers

from ..utility.io import search_sbo_label
from ..utility.entities import reaction_equation_to_dict

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# investigate with memote
# -----------------------

def run_memote(model: cobra.Model, type:Literal['json','html']='html', 
               return_res:bool=False, save_res:str|None=None, verbose:bool=False) -> dict|str|None:
    """Run the memote snapshot function on a given model loaded with COBRApy.

    Args:
        - model (cobra.Model): 
            The model loaded with COBRApy.
        - type (Literal['json','html'], optional): 
            Type of report to produce. 
            Can be 'html' or 'json'. 
            Defaults to 'html'.
        - return_res (bool, optional): 
            Option to return the result. 
            Defaults to False.
        - save_res (str | None, optional): 
            If given a path string, saves the report
            under the given path. Defaults to None.
        - verbose (bool, optional): 
            Produce a more verbose ouput. 
            Defaults to False.

    Raises:
        ValueError: Unknown input for parameter type

    Returns:
        dict|str|None: The json dictionary, the html string or none.
    """

    # verbose output I
    if verbose:
        print('\n# -------------------\n# Analyse with MEMOTE\n# -------------------')
        start = time.time()

    # run memote
    ret, res = memote.suite.api.test_model(model, sbml_version=None, results=True,
                                           pytest_args=None, exclusive=None, skip=None, 
                                           experimental=None, solver_timeout=10)
    
    # load depending on type 
    match type:
        case 'html':
            snap = memote.suite.api.snapshot_report(res, html=True)
            result = snap
        case 'json':
            snap = memote.suite.api.snapshot_report(res, html=False)
            result = json.loads(snap)
        case _:
            message = f'Unknown input for parameter type: {type} '
            raise ValueError(message)
        
    # option to save report
    if save_res:
        with open(save_res, 'w') as f:
            f.write(result)

    # verbose output II
    if verbose:
        end = time.time()
        print(F'\ttotal time: {end - start}s')

    # option to return report
    if return_res:
        return result
    

def get_memote_score(memote_report: dict) -> float:
    """Extracts MEMOTE score from report

    Args:
        - memote_report (dict): Output from run_memote.

    Returns:
        float: MEMOTE score
    """
    return memote_report['score']['total_score']


# get basic model info
# --------------------

def get_num_reac_with_gpr(model:cobra.Model) -> int:
    """Extract the number of reactions that have a gene production rule
    from a given model.

    Args:
        model (cobra.Model): The model loaded with COBRApy.

    Returns:
        int: The number of reactions with a GPR.
    """

    reac_with_gpr = 0
    for reac in model.reactions:
        # check for GPR
        if len(reac.genes) > 0:
            reac_with_gpr += 1

    return reac_with_gpr


def get_orphans_deadends_disconnected(model: cobraModel) -> tuple[list[str], list[str], list[str]]:
    """Uses MEMOTE functions to extract orphans, deadends and disconnected metabolites

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        tuple: Lists of metabolites that might cause errors (1) - (3) 
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


# @TODO: what about exchange reactions not starting with an EX - as usually the case within the BiGG namespace? 
def get_mass_charge_unbalanced(model: cobraModel) -> tuple[list[str], list[str]]:
    """Creates lists of mass and charge unbalanced reactions,vwithout exchange reactions since they are unbalanced per definition

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        tuple: Lists of reactions that might cause errors (1) & (2)
        (1) list: List of mass unbalanced reactions
        (2) list: List of charge unbalanced reactions
    """

    mass_unbalanced = consistency.find_mass_unbalanced_reactions(
        model.reactions)
    charge_unbalanced = consistency.find_charge_unbalanced_reactions(
        model.reactions)

    mass_list = []
    if len(mass_unbalanced) > 0:
        for reac in mass_unbalanced:
            if (reac.id[:2] != 'EX'):
                mass_list.append(reac.id)

    charge_list = []
    if len(charge_unbalanced) > 0:
        for reac in charge_unbalanced:
            if (reac.id[:2] != 'EX'):
                charge_list.append(reac.id)

    return mass_list, charge_list


# other
# -----

def get_metabs_with_one_cvterm(model: libModel) -> list[str]:
    """Reports metabolites which have only one annotation, can be used as basis for further annotation research

    Args:
        - model (libModel): Model loaded with libSBML

    Returns:
        list: Metabolite Ids with only one annotation
    """
    spe = model.getListOfSpecies()

    only_one = [] #safe metab with only BiGG annotation
    for sb in spe:
        if sb.isSetId():
            pid = sb.getId()
            if sb.getCVTerm(0).getNumResources() == 1:
                only_one.append(pid)
                
    return only_one

def get_reactions_per_sbo(model: libModel) -> dict:
    """Counts number of reactions of all SBO Terms present

    Args:
        - model (libModel): Model loaded with libSBML

    Returns:
        dict: SBO Term as keys and number of reactions as values
    """
    sbos_dict = {}
    for react in model.getListOfReactions():
        sbo = react.getSBOTerm()
        if sbo in sbos_dict.keys():
            sbos_dict[sbo] += 1
        else: 
            sbos_dict[sbo] = 1
    return sbos_dict

def plot_rea_sbo_single(model: libModel):
    """Plots reactions per SBO Term in horizontal bar chart

    Args:
        - model (libModel): Model loaded with libSBML

    Returns:
        plot: Pandas Barchart
    """
    df = pd.DataFrame(get_reactions_per_sbo(model), index=[0]).T.reset_index().rename({0:model.id, 'index': 'SBO-Term'}, axis=1)
    df = df[df[model.id]>3]
    df['SBO-Name'] = df['SBO-Term'].apply(search_sbo_label)
    fig = df.drop('SBO-Term', axis=1).sort_values(model.id).set_index('SBO-Name').plot.barh(width=.8, figsize=(8,10))
    fig.set_ylabel('')
    fig.set_xlabel('number of reactions', fontsize=16)
    fig.legend(loc='lower right')
    return fig