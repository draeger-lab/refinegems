#!/usr/bin/env python
""" Provides functions to simulate growth on any medium

Tailored to work with the Synthetic Nasal Medium (SNM), should
work with any medium as long as its defined in a csv with ; as 
delimiter and BiGG Ids for the compounds.
"""

import cobra
import pandas as pd
import numpy as np
from libsbml import Model

__author__ = "Famke Baeuerle"

def get_default_uptake(model):
    """Determines which metabolites are used in the standard medium

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        list: metabolites consumed in standard medium
    """
    fluxes = model.optimize().fluxes
    default_uptake = []
    for index, value in fluxes.items():
        if "EX" in index:
            if value < 0:
                default_uptake.append(index)
    return default_uptake

def load_medium_custom(mediumpath):
    """Helper function to read medium csv

    Args:
        mediumpath (Str): path to csv file with medium

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium['BiGG_R']='R_EX_'+ medium['BiGG']+ '_e'
    medium['BiGG_EX']='EX_'+ medium['BiGG']+ '_e'
    return medium

def load_medium_from_db(mediumpath, mediumname):
    """Helper function to read standard media_db.csv

    Args:
        mediumpath (Str): path to csv file with medium database
        mediumname (Str): name of medium to test growth on

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium = medium.loc[medium['medium'] == mediumname]
    medium['BiGG_R']='R_EX_'+ medium['BiGG']+ '_e'
    medium['BiGG_EX']='EX_'+ medium['BiGG']+ '_e'
    return medium
    
def find_missing_exchanges(model, medium):
    """Look for exchange reactions needed by the medium but not in the model

    Args:
        model (libsbml-model): model loaded with libsbml
        medium (df): dataframe of medium csv

    Returns:
        list: all exchanges missing in the model but given in medium
    """
    medium_list = medium['BiGG_R'].dropna().tolist()
    missing_exchanges = []
    for exchange in medium_list:
        if model.getReaction(exchange) == None:
            missing_exchanges.append(exchange)
    return missing_exchanges

def modify_medium(medium, missing_exchanges):
    """Remove exchanges from medium that are not in the model

    Args:
        medium (df): dataframe of medium csv
        missing_exchanges (list): exchanges not in the model

    Returns:
        dict: medium that can be used with the model
    """
    medium_dict = dict.fromkeys(medium['BiGG_EX'].dropna().tolist(), 10.0)
    for exchange in missing_exchanges:
        exchange = exchange[2:]
        medium_dict.pop(exchange)
    return medium_dict

def find_missing_essential(model, medium, default_uptake):
    """Report which exchange reactions are needed for growth
        combines default uptake and valid new medium

    Args:
        model (cobra-model): model loaded with cobrapy
        medium (dict): medium with exchanges that are present in the model
        default_uptake (list): metabolites consumed in standard medium

    Returns:
        list: all metabolites which lead to zero growth if blocked
    """
    medium.update({ i : 10.0 for i in default_uptake })
    med = model.medium
    model.medium = medium
    essential = []
    for metab in medium.keys():
        model.reactions.get_by_id(metab).lower_bound = 0
        sol = model.optimize()
        if sol.objective_value < 1e-9: # and sol.objective_value > -1e-9: # == 0 no negative growth!
            essential.append(metab)
        model.reactions.get_by_id(metab).lower_bound = -10
    return essential

def find_minimum_essential(medium, essential):
    """Report metabolites necessary for growth and not in custom medium

    Args:
        medium (df): dataframe of medium csv
        essential (list): metabolites which lead to zero growth if blocked

    Returns:
        list: metabolites not present in the medium but necessary for growth
    """
    minimum = []
    for metab in essential:
        if metab not in medium['BiGG_EX'].tolist():
            minimum.append(metab)
    return minimum

def simulate_minimum_essential(model, medium, minimum):
    """Simulate growth with custom medium plus necessary uptakes

    Args:
        model (cobra-model): model loaded with cobrapy
        medium (dict): medium with exchanges present in the model
        minimum (list): metabolites not present in the medium but necessary for growth

    Returns:
        float: growth value
    """
    medium.update({ i : 10.0 for i in minimum })
    try:
        if (medium['EX_o2_e'] == 10.0):
            medium['EX_o2_e'] = 20.0
    except KeyError:
        pass
    med = model.medium
    model.medium = medium
    sol = model.optimize()
    return sol.objective_value

def growth_simulation(model_cobra, model_libsbml, mediumpath, mediumname):
    """Executes all steps to determine growth on custom medium

    Args:
        model_cobra (cobra-model): model loaded with cobrapy
        model_libsbml (libsbml-model): model loaded with libsbml
        mediumpath (Str): path to csv file with medium

    Returns:
        tuple: (added exchanges, exchanges missing in model, growth value)
    """
    medium = load_medium_from_db(mediumpath, mediumname)
    default_uptake = get_default_uptake(model_cobra)
    missing_exchanges = find_missing_exchanges(model_libsbml, medium)
    medium_dict = modify_medium(medium, missing_exchanges)
    essential = find_missing_essential(model_cobra, medium_dict, default_uptake)
    minimum = find_minimum_essential(medium, essential)
    # reload medium definition to exclude all previous additions
    medium_dict = modify_medium(medium, missing_exchanges)
    growth_value = simulate_minimum_essential(model_cobra, medium_dict, minimum)
    doubling_time = (np.log(2)/growth_value) * 60
    return minimum, missing_exchanges, growth_value, doubling_time