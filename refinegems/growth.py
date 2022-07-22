#!/usr/bin/env python
""" Provides functions to simulate growth on any medium

Tailored to work with the Synthetic Nasal Medium (SNM), should
work with any medium as long as its defined in a csv with ; as
delimiter and BiGG Ids for the compounds.

Outputs a table where

essential = metabolites not present in the medium but necessary for growth

missing = all exchanges missing in the model but given in medium
"""

import pandas as pd
import numpy as np
from refinegems.load import load_all_media_from_db
from cobra.medium import minimal_medium

__author__ = "Famke Baeuerle"


def get_default_uptake(model):
    """Determines which metabolites are used in the standard medium

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        list: metabolites consumed in standard medium
    """
    with model:
        sol = model.optimize()
        fluxes = sol.fluxes
        default_uptake = []
        for index, value in fluxes.items():
            if "EX" in index:
                if value < 0:
                    default_uptake.append(index)
    return default_uptake


def get_minimum_default(model):
    with model:
        #new_medium = {i: 10.0 for i in default_uptake}
        #model.medium = new_medium
        minimal = minimal_medium(model)
        uptake = []
        for index, value in minimal.items():
            if value < 0:
                uptake.append(index)
    return list(minimal.index)


def get_default_secretion(model):
    """Checks fluxes after fba, if positive the metabolite is produced

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        list: BiGG Ids of produced metabolites
    """
    with model:
        fluxes = model.optimize().fluxes
        default_secretion = []
        for index, value in fluxes.items():
            if value > 0:
                default_secretion.append(index)
    return default_secretion


def get_missing_exchanges(model, medium):
    """Look for exchange reactions needed by the medium but not in the model

    Args:
        model (cobra-model): model loaded with cobrapy
        medium (df): dataframe of medium csv

    Returns:
        list: all exchanges missing in the model but given in medium
    """
    medium_list = medium['BiGG_EX'].dropna().tolist()
    missing_exchanges = []
    for exchange in medium_list:
        try:
            model.reactions.get_by_id(exchange)
        except(KeyError):
            missing_exchanges.append(exchange)
    return missing_exchanges


def modify_medium(medium, missing_exchanges):
    """Remove exchanges from medium that are not in the model

    Args:
        medium (df): dataframe of medium csv
        missing_exchanges (list): exchanges not in the model

    Returns:
        dict: growth medium that can be used with the model
    """
    # needed or else KeyError later
    medium_dict = dict.fromkeys(medium['BiGG_EX'].dropna().tolist(), 10.0)
    for exchange in missing_exchanges:
        medium_dict.pop(exchange)
    return medium_dict


def find_missing_essential(model, medium, default_uptake):
    """Report which exchange reactions are needed for growth
        combines default uptake and valid new medium

    Args:
        model (cobra-model): model loaded with cobrapy
        medium (dict): custom medium with exchanges that are present in the model
        default_uptake (list): metabolites consumed in standard medium

    Returns:
        list: all metabolites which lead to zero growth if blocked
    """
    with model:
        default_medium = {i: 10.0 for i in default_uptake}
        new_medium = {**medium, **default_medium}
        model.medium = new_medium
        essential = []
        for metab in new_medium.keys():
            with model:
                model.reactions.get_by_id(metab).lower_bound = 0
                sol = model.optimize()
                if sol.objective_value < 1e-5:  # and sol.objective_value > -1e-9: # == 0 no negative growth!
                    essential.append(metab)
                model.reactions.get_by_id(metab).lower_bound = -10
    return essential


def find_minimum_essential(medium, essential):
    """Report metabolites necessary for growth and not in custom medium

    Args:
        medium (df): dataframe of medium csv
        essential (list): echanges of metabolites which lead to zero growth if blocked

    Returns:
        list: metabolites not present in the medium but necessary for growth
    """
    minimum = []
    for metab in essential:
        if metab not in medium['BiGG_EX'].tolist():
            minimum.append(metab)
    return minimum


def simulate_minimum_essential(model, growth_medium, minimum):
    """Simulate growth with custom medium plus necessary uptakes

    Args:
        model (cobra-model): model loaded with cobrapy
        growth_medium (dict): growth medium with exchanges present in the model
        minimum (list): exchanges of metabolites not present in the medium but necessary for growth

    Returns:
        float: growth value
    """
    with model:
        min_medium = {i: 10.0 for i in minimum}
        new_medium = {**growth_medium, **min_medium}
        try:
            if (new_medium['EX_o2_e'] == 10.0):
                new_medium['EX_o2_e'] = 20.0
        except KeyError:
            print('No Oxygen Exchange Reaction')
            pass
        model.medium = new_medium
        sol = model.optimize()
    return sol.objective_value


def get_all_minimum_essential(model, mediumpath):
    """Returns metabolites necessary for growth and not in media

    Args:
        model (cobra-model): model loaded with cobrapy
        mediumpath (string): path to csv with media definitions

    Returns:
        DataFrame: information on different media which metabs are missing
    """
    media = load_all_media_from_db(mediumpath)
    default_uptake = get_default_uptake(model)
    mins = pd.DataFrame()
    for medium in media:
        missing_exchanges = get_missing_exchanges(model, medium)
        medium_dict = modify_medium(medium, missing_exchanges)
        essential = find_missing_essential(model, medium_dict, default_uptake)
        minimum = find_minimum_essential(medium, essential)
        mins[medium['medium'][0]] = pd.Series(minimum)
    return mins


def get_growth_one_medium(model, medium):
    """Simulates growth on given medium

    Args:
        model (cobra-model): model loaded with cobrapy
        medium (pandas-DataFrame): table containing metabolites present in the medium

    Returns:
        DataFrame: information on growth behaviour on medium
    """
    default_uptake = get_default_uptake(model)
    minimal_medium = get_minimum_default(model) # use this instead of default_uptake
    missing_exchanges = get_missing_exchanges(model, medium)
    medium_dict = modify_medium(medium, missing_exchanges)
    essential = find_missing_essential(model, medium_dict, minimal_medium)#default_uptake)
    minimum = find_minimum_essential(medium, essential)

    medium_dict = modify_medium(medium, missing_exchanges)
    growth_value = simulate_minimum_essential(model, medium_dict, minimum)
    doubling_time = (np.log(2) / growth_value) * 60
    exchanges = [[medium['medium'][0]], minimum,
                 missing_exchanges, [growth_value], [doubling_time]]
    df_growth = pd.DataFrame(exchanges,
                             ['medium',
                              'essential',
                              'missing',
                              'growth_value',
                              'doubling_time [min]']).T
    return df_growth


def get_growth_selected_media(model, mediumpath, media):
    """Simulates growth on all selected media

    Args:
        model (cobra-model): model loaded with cobrapy
        mediumpath (string): path to csv with media definitions
        media (list): media to simulate on (must be in csv)

    Returns:
        DataFrame: information on growth behaviour on selected media
    """
    all_media = load_all_media_from_db(mediumpath)
    selected_media = [x for x in all_media if x['medium'][0] in media]
    growth = pd.DataFrame()
    for medium in selected_media:
        growth = growth.append(
            get_growth_one_medium(model, medium), 
            ignore_index=True)
    return growth


def get_essential_reactions(model):
    """Knocks out each reaction, if no growth is detected the reaction is seen as essential

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        list: BiGG Ids of essential reactions
    """
    ess = []
    for reaction in model.reactions:
        with model as model:
            reaction.knock_out()
            model.optimize()
            if model.objective.value <= 11:
                print('%s blocked (bounds: %s), new growth rate %f $' %
                    (reaction.id, str(reaction.bounds), model.objective.value))
                ess.append(reaction.id)
                
    #### other implementation via bounds ###
    # medium = model.medium
    # ess = []
    # for content in medium.keys():
    #     model.reactions.get_by_id(content).lower_bound = 0.0
    #     solution = model.optimize().objective_value
    #     if solution < 1e-9:
    #         ess.append(content)
    #     model.reactions.get_by_id(content).lower_bound = -10.0

    return ess