#!/usr/bin/env python
""" Provides functions to investigate the model and test with memote

These functions enable simple testing of any model using memote and 
access to its number of reactions, metabolites and genes.
"""

import os
import memote
import json
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
