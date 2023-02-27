#!/usr/bin/env python
""" Provides functions to investigate the model and test with memote

These functions enable simple testing of any model using memote and
access to its number of reactions, metabolites and genes.
"""

import os
import memote
import json
import pandas as pd
import numpy as np
from cobra import Reaction
from memote.support import consistency
# needed by memote.support.consitency
from memote.support import consistency_helpers as con_helpers
from refinegems.io import load_model_cobra, load_model_libsbml, search_sbo_label

__author__ = "Famke Baeuerle"

DISSIPATION_RXNS = {
    'ATP':'atp_c + h2o_c --> adp_c + h_c + pi_c',
    'CTP':'ctp_c + h2o_c --> cdp_c + h_c + pi_c',
    'GTP':'gtp_c + h2o_c --> gdp_c + h_c + pi_c',
    'UTP':'utp_c + h2o_c --> udp_c + h_c + pi_c',
    'ITP':'itp_c + h2o_c --> idp_c + h_c + pi_c',
    'NADH':'nadh_c --> h_c + nad_c',
    'NADPH':'nadph_c --> h_c + nadp_c',
    'FADH2':'fadh2_c --> 2 h_c + fad_c',
    'FMNH2':'fmnh2_c --> 2 h_c + fmn_c',
    'Q8H2':'q8h2_c --> 2 h_c + q8_c',
    'MQL8':'mql8_c --> 2 h_c + mqn8_c',
    'DMMQL8':'2dmmql8_c --> 2 h_c + 2dmmq8_c',
    'ACCOA':'h2o_c + accoa_c --> h_c + ac_c + coa_c',
    'GLU':'glu__L_c + h2o_c --> akg_c + nh4_c + 2 h_c',
    'PROTON':'h_p --> h_c'
    }


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
                                           pytest_args=None, exclusive=None, skip=None, 
                                           experimental=None, solver_timeout=10)
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
        tuple: (list of orphans, deadends, disconnected metabolites)
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


def get_model_info(modelpath):
    """Reports core information of given model

    Args:
        modelpath (string): path to model file

    Returns:
        DataFrame: overview on model parameters
    """
    model_libsbml = load_model_libsbml(modelpath)
    model_cobra = load_model_cobra(modelpath)
    name, reac, metab, genes = initial_analysis(model_libsbml)
    orphans, deadends, disconnected = get_orphans_deadends_disconnected(
        model_cobra)
    mass_unbal, charge_unbal = get_mass_charge_unbalanced(model_cobra)
    model_info = pd.DataFrame([name,
                               reac,
                               metab,
                               genes,
                               orphans,
                               deadends,
                               disconnected,
                               mass_unbal,
                               charge_unbal],
                              ['model',
                               '#reactions',
                               '#metabolites',
                               '#genes',
                               'orphans',
                               'deadends',
                               'disconnected',
                               'mass unbalanced',
                               'charge unbalanced']).T

    return model_info

def parse_reaction(eq, model): # from alina
    """Parses a reaction equation string to dictionary

    Args:
        eq (string): Equation of a reaction
        model (cobra-model): model loaded with cobrapy

    Returns:
       dict: metabolite ids as keys and their coefficients as values (negative = educts, positive = products)
    """
    eq = eq.split(' ')
    eq_matrix={}
    are_products = False
    coeff = 1
    for i,part in enumerate(eq):
        if part == '-->':
            are_products = True
            continue          
        if part == '+':
            continue
        if part == '2':
            coeff = 2
            continue
        if are_products:
            eq_matrix[model.metabolites.get_by_id(part)] = 1*coeff
            coeff = 1
        else:
            eq_matrix[model.metabolites.get_by_id(part)] = -1*coeff
            coeff = 1
    return eq_matrix

def get_egc(model):
    """Energy-generating cycles represent thermodynamically infeasible states. Charging of energy metabolites without any energy source causes such cycles. Detection method is based on (Fritzemeier et al., 2017)

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        df: table with possible EGCs
    """
    dissipation_rxns = pd.DataFrame(DISSIPATION_RXNS.items(), columns=['type','equation'])
    
    with model: 
    # add dissipation reactions
        for i, row in dissipation_rxns.iterrows():
            try:
                met_atp = parse_reaction(row['equation'], model)
                rxn = Reaction(row['type'])
                rxn.name = 'Test ' + row['type'] + ' dissipation reaction'
                rxn.add_metabolites(met_atp)
                model.add_reaction(rxn)
            except(KeyError):
                dissipation_rxns.drop(dissipation_rxns[dissipation_rxns['type'] == row['type']].index, inplace=True)
            
        for rxn in model.reactions:
            if 'EX_' in rxn.id:
                rxn.lower_bound = 0.0
                rxn.upper_bound = 0.0
                #print('Set exchange rxn to 0', rxn.name)
            # set reversible reactions fluxes to [-1,1]    
            elif rxn.reversibility: 
                rxn.lower_bound = -1.0
                rxn.upper_bound = 1.0
                #print('Reversible rxn', rxn.name)
            # irreversible reactions have fluxes [0.1]    
            else:
                rxn.lower_bound = 0.0
                rxn.upper_bound = 1.0
                #print('Irreversible rxn', rxn.name)
                
        df_fluxes = pd.DataFrame()
        objval = dict()
        # optimize by choosing one of dissipation reactions as an objective
        for i, row in dissipation_rxns.iterrows():
            model.objective = row['type']
            #print('Set objective to', row['type'], ':', model.optimize().objective_value)
            objval[row['type']] = model.optimize().objective_value
            if model.optimize().objective_value > 0.0:
                df = pd.DataFrame.from_dict([model.optimize().fluxes]).T.replace(0, np.nan).dropna(axis=0).rename({'fluxes':row['type']}, axis=1)
                df_fluxes = pd.concat([df_fluxes, df], axis=1)
            else:
                df_fluxes[row['type']] = np.nan
        df_fluxes = pd.concat([df_fluxes,pd.DataFrame.from_dict([objval])])
    return df_fluxes.T.reset_index().rename({'index':'BOF', 0:'objective value'}, axis=1).fillna('')

def get_metabs_with_one_cvterm(model):
    """reports metabolites which have only one annotation, 
    can be used as basis for further annotation research

    Args:
        model (libsbml-model): model loaded with libsbml

    Returns:
        list: metabolite Ids
    """
    spe = model.getListOfSpecies()

    only_one = [] #safe metab with only BiGG annotation
    for sb in spe:
        if sb.isSetId():
            pid = sb.getId()
            if sb.getCVTerm(0).getNumResources() == 1:
                only_one.append(pid)
                
    return only_one

def get_reactions_per_sbo(model):
    sbos_dict = {}
    for react in model.getListOfReactions():
        sbo = react.getSBOTerm()
        if sbo in sbos_dict.keys():
            sbos_dict[sbo] += 1
        else: 
            sbos_dict[sbo] = 1
    return sbos_dict

def get_sbo_plot_single(model):
    df = pd.DataFrame(get_reactions_per_sbo(model), index=[0]).T.reset_index().rename({0:model.id, 'index': 'SBO-Term'}, axis=1)
    df = df[df[model.id]>3]
    df['SBO-Name'] = df['SBO-Term'].apply(search_sbo_label)
    fig = df.drop('SBO-Term', axis=1).sort_values(model.id).set_index('SBO-Name').plot.barh(width=.8, figsize=(8,10))
    fig.set_ylabel('')
    fig.set_xlabel('number of reactions', fontsize=16)
    fig.legend(loc='lower right')
    return fig