#!/usr/bin/env python
""" Provides functions to investigate the model and test with memote

These functions enable simple testing of any model using memote and 
access to its number of reactions, metabolites and genes.
"""

import os
import memote
import json
from memote.support import consistency
# needed by memote.support.consitency
from memote.support import consistency_helpers as con_helpers 

__author__ = "Famke Baeuerle"


def run_memote_sys(modelfile):
    """run memote on linux machine

    Args:
        modelfile (cobra-model): model loaded with cobrapy
    """
    cmd = 'memote report snapshot ' + str(modelfile)
    os.system(cmd)


def run_memote(model):
    """runs memote to obtain total score

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        float: memote score of model
    """
    ret, res = memote.suite.api.test_model(model, sbml_version=None, results=True,
                                            pytest_args=None, exclusive=None, skip=None, experimental=None, solver_timeout=10)
    snap = memote.suite.api.snapshot_report(res, html=False)
    result = json.loads(snap)
    totalScore = result['score']['total_score']
    return totalScore



def initial_analysis(model):
    """extracts most important numbers of GEM

    Args:
        model (libsbml-model): model loaded with libsbml

    Returns:
        tuple: (name of model, no of reactions, no of metabolites, no of genes)
    """
    name = model.getId()
    reactions = model.getNumReactions()
    metabolites = model.getNumSpecies()
    genes = len(model.getPlugin(0).getListOfGeneProducts())
    return name, reactions, metabolites, genes


def get_orphans_deadends_disconnected(model):
    """Uses memote functions to extract orphans, deadends and disconnected metabolites

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        tuple: (list of BiGG ids of orphans, deadends, disconnected metabolites)
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

def get_mass_charge_unbalanced(model):
    """creates lists of mass and charge unbalanced reactions, 
       without exchange reactions since they are unbalanced per definition

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        tuple: (list of mass unbalanced, charge unbalanced reactions)
    """
    
    mass_unbalanced = consistency.find_mass_unbalanced_reactions(model.reactions)
    charge_unbalanced = consistency.find_charge_unbalanced_reactions(model.reactions)
    
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
