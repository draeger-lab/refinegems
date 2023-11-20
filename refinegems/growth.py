#!/usr/bin/env python
""" Provides functions to simulate growth on any medium

Tailored to work with the media denoted in the local db, should work with any medium as long as its defined in a csv with ; as delimiter and BiGG Ids for the compounds. Use refinegems.io.load_medium_custom and hand this to the growth_one_medium_from_default or growth_one_medium_from_minimum function.
"""

import logging
import pandas as pd
import numpy as np
from refinegems.io import load_medium_from_db_for_growth
from cobra import Reaction
from cobra import Model as cobraModel
import cobra
import re

__author__ = "Famke Baeuerle and Carolin Brune"



# @TEST
def set_bounds_to_default(model: cobraModel, reac_bounds = None):
    """Set the reactions bounds of a model to given default values.
    (Ir)reversibility is retained.

    Args:
        model (cobraModel): The model loaded with COBRApy.
        reac_bounds ([NoneType, str, tuple[float]], optional): The setting for the new reaction bounds. 
            Defaults to None. If None or "cobra", uses the COBRApy in-built default values (-1000.0, 1000.0).
            The user can set personal values by entering a tuple of two floats.

    Raises:
        ValueError: Problematic input for bounds, if neither None, "cobra" or a tuple of floats is entered for reac_bounds.
    """

    # user-specific default bounds (tuple of two floats)
    if (type(reac_bounds) is tuple 
         and len(reac_bounds) == 2 
         and type(reac_bounds[0]) is float 
         and type(reac_bounds[1]) is float):
        pass
    # use COBRApy-internal default bounds (None or string 'cobra')
    elif reac_bounds is None or (type(reac_bounds) is str and reac_bounds == 'cobra'):
        reac_bounds = cobra.Configuration().bounds
    # unknown input
    else:
        raise ValueError(f'Problematic input for bounds: {reac_bounds}')
    
    # apply the bounds to the model
    for reaction in model.reactions:

        # reactions is originally disabled
        if reaction.upper_bound == 0.0 and reaction.lower_bound == 0.0:
            pass
        # forward only
        elif reaction.lower_bound == 0.0:
            reaction.upper_bound = reac_bounds[1]
        # backward only
        elif reaction.upper_bound == 0.0:
            reaction.lower_bound = reac_bounds[0]
        # reversible or broken
        else:
            reaction.bounds = reac_bounds


# @TEST
def get_uptake(model: cobraModel, type: str, exchange_regex='^EX') -> list[str]:
    """Compute the list of exchange reactions that have fluxes > 0 under certain conditions.

    Args:
        model (cobraModel): A cobra Model to be tested.
        type (str): Type of uptake, can be 'minimal'/'min' or 'standard'/'std'.
        exchange_regex (str, optional): Regex-compatible string to determine exchange reactions. Defaults to '^EX'.

    Raises:
        ValueError: Unknown type for uptake, if type not in ['minimal','min','standard','std']

    Returns:
        list[str]: List of non-zero flux exchange reactions under the set type.
    """

    match type:
        # return minimal 
        case 'minimal' | 'min':
            with model:
                minimal = cobra.medium.minimal_medium(model)
                print(minimal)
                return list(minimal.index)
        # return standart, non-zero flux compounds
        case 'standard' | 'std':
            with model:
                sol = model.optimize()
                fluxes = sol.fluxes
                uptake = []
                regexp = re.compile(exchange_regex)
                for index, value in fluxes.items():
                    if regexp.search(index):
                        if value < 0:
                            uptake.append(index)
            return uptake
        case _:
            raise ValueError(f'Unknown type for uptake: {type}')



def get_secretion(model: cobraModel) -> list[str]:
    """Returns the list of exchange reactions for compounds that are secreted in the current version of the model.

    Args:
        model (cobraModel): The cobra model to be tested.

    Returns:
        list[str]: The list of IDs of secretion reactions
    """

    with model:

        sf = model.summary().secretion_flux
        s = sf[sf['flux'] < 0.0].index.tolist()

    return s
    
def get_production(model: cobraModel) -> list[str]:
    """Checks fluxes after FBA, if positive the metabolite is produced.

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: Ids of produced metabolites
    """

    with model:
        fluxes = model.optimize().fluxes
        production = []
        for index, value in fluxes.items():
            if value > 0:
                production.append(index)
    return production


# combine with below
def get_missing_exchanges(model: cobraModel, medium: pd.DataFrame) -> list[str]:
    """Look for exchange reactions needed by the medium but not in the model

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - medium (pd.DataFrame): Dataframe with medium definition

    Returns:
        list[str]: Ids of all exchanges missing in the model but given in medium
    """
    medium_list = medium['BiGG_EX'].dropna().tolist()
    missing_exchanges = []
    for exchange in medium_list:
        try:
            model.reactions.get_by_id(exchange)
        except(KeyError):
            missing_exchanges.append(exchange)
    return missing_exchanges


# already done with medium.py add_medium_to_model
def modify_medium(medium: pd.DataFrame, missing_exchanges: list[str]) -> dict:
    """Helper function: Remove exchanges from medium that are not in the model to avoid KeyError

    Args:
        - medium (pd.DataFrame): Dataframe with medium definition
        - missing_exchanges (list): Ids of exchanges not in the model

    Returns:
        dict: Growth medium definition that can be used with the model (f.ex {'EX_glc__D_e' : 10.0})
    """
    growth_medium = dict.fromkeys(medium['BiGG_EX'].dropna().tolist(), 10.0)
    for exchange in missing_exchanges:
        growth_medium.pop(exchange)
    return growth_medium


# idea: remove default uptake from params list
# recheck connection to add_medium_to_model from medium.py
# where to elegantly put anaerobic option
# not exhaustive, might lead to problems
# rename?
def find_missing_essential(model: cobraModel, growth_medium: dict, default_uptake: list[str], anaerobic: bool) -> list[str]:
    """Report which exchange reactions are needed for growth, combines default uptake and valid new medium

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - growth_medium (dict): Growth medium definition that can be used with the model. Output of modify_medium.
        - default_uptake (list[str]): Metabolites consumed in standard medium
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        list[str]: Ids of exchanges of all metabolites which lead to zero growth if blocked
    """
    with model:
        default_medium = {i: 10.0 for i in default_uptake}
        if anaerobic and ('EX_o2_e' in default_medium): default_medium['EX_o2_e'] = 0.0
        if anaerobic and ('EX_o2_e' in growth_medium): growth_medium['EX_o2_e'] = 0.0
        new_medium = {**growth_medium, **default_medium}
        try:
            model.medium = new_medium
        except(ValueError):
            logging.info('Change upper bounds to 1000.0 and lower bounds to -1000.0 to make model simulatable.')
            # ..................................
            #for reaction in model.reactions:
            #    set_bounds_to_default(reaction)
            # replace with new
            set_bounds_to_default(model)
            # ..................................
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

# combine with above, only occur together
def find_minimum_essential(medium: pd.DataFrame, essential: list[str]) -> list[str]:
    """Report metabolites necessary for growth and not in custom medium

    Args:
        - medium (pd.DataFrame): Dataframe with medium definition
        - essential (list[str]): Ids of all metabolites which lead to zero growth if blocked. Output of find_missing_essential.

    Returns:
        list[str]: Ids of exchanges of metabolites not present in the medium but necessary for growth
    """
    minimum = []
    for metab in essential:
        if metab not in medium['BiGG_EX'].tolist():
            minimum.append(metab)
    return minimum

# anaerobic see above
# rename: simulate_growth_with essential_supplements()
# verallgemeinern?
def simulate_minimum_essential(model: cobraModel, growth_medium: dict, minimum: list[str], anaerobic: bool) -> float:
    """Simulate growth with custom medium plus necessary uptakes

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - growth_medium (dict): Growth medium definition that can be used with the model. Output of modify_medium.
        - minimum (list[str]): Ids of exchanges of metabolites not present in the medium but necessary for growth. Output of find_minimum_essential.
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        float: Growth value in mmol per (gram dry weight) per hour
    """
    with model:
        min_medium = {i: 10.0 for i in minimum}
        new_medium = {**growth_medium, **min_medium}
        try:
            if (new_medium['EX_o2_e'] == 10.0):
                new_medium['EX_o2_e'] = 20.0 if not anaerobic else 0.0
        except KeyError:
            print('No Oxygen Exchange Reaction')
            pass
        try:
            model.medium = new_medium
        except(ValueError):
            logging.info('Change upper bounds to 1000.0 and lower bounds to -1000.0 to make model simulatable.')
            # ..................................
            # changed old to new
            # for reaction in model.reactions:
            #     set_bounds_to_default(reaction)
            
            set_bounds_to_default(model)
            # ..................................
            model.medium = new_medium
        sol = model.optimize()
    return sol.objective_value

# entkoppelt
# ??????????
def get_all_minimum_essential(model: cobraModel, media: list[str]) -> pd.DataFrame:
    """Returns metabolites necessary for growth and not in media

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - media (list[str]): Containing the names of all media for which the growth essential metabolites not contained in the media should be returned

    Returns:
        pd.DataFrame: information on different media which metabs are missing
    """
    default_uptake = get_uptake(model,'std')
    mins = pd.DataFrame()
    for medium in media:
        medium_df = load_medium_from_db_for_growth(medium)
        missing_exchanges = get_missing_exchanges(model, medium_df)
        medium_dict = modify_medium(medium_df, missing_exchanges)
        essential = find_missing_essential(model, medium_dict, default_uptake)
        minimum = find_minimum_essential(medium_df, essential)
        mins[medium['medium'][0]] = pd.Series(minimum)
    return mins

# combine with below
# growth_one_medium
# output?
def growth_one_medium_from_default(model: cobraModel, medium: pd.DataFrame, anaerobic: bool) -> pd.DataFrame:
    """Simulates growth on given medium, adding missing metabolites from the default uptake

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - medium (pd.DataFrame): Dataframe with medium definition
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        pd.DataFrame: Information on growth behaviour on given medium
    """
    default_uptake = get_uptake(model,'std')
    missing_exchanges = get_missing_exchanges(model, medium)
    medium_dict = modify_medium(medium, missing_exchanges)
    essential = find_missing_essential(model, medium_dict, default_uptake, anaerobic)
    minimum = find_minimum_essential(medium, essential)

    medium_dict = modify_medium(medium, missing_exchanges)
    growth_value = simulate_minimum_essential(model, medium_dict, minimum, anaerobic)
    doubling_time = (np.log(2) / growth_value) * 60
    medium_name = medium['medium'][0] if not anaerobic else f'{medium["medium"][0]}[-O2]'
    exchanges = [[medium_name], minimum,
                 missing_exchanges, [growth_value], [doubling_time]]
    df_growth = pd.DataFrame(exchanges,
                             ['medium',
                              'essential',
                              'missing exchanges',
                              'growth_value',
                              'doubling_time [min]']).T
    return df_growth


def growth_one_medium_from_minimal(model: cobraModel, medium: pd.DataFrame, anaerobic: bool) -> pd.DataFrame:
    """Simulates growth on given medium, adding missing metabolites from a minimal uptake

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - medium (pd.DataFrame): Dataframe with medium definition
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        pd.DataFrame: Information on growth behaviour on given medium
    """
    minimal_uptake = get_uptake(
        model, 'min')  # use this instead of default_uptake
    missing_exchanges = get_missing_exchanges(model, medium)
    medium_dict = modify_medium(medium, missing_exchanges)
    essential = find_missing_essential(model, medium_dict, minimal_uptake, anaerobic)
    minimum = find_minimum_essential(medium, essential)

    medium_dict = modify_medium(medium, missing_exchanges)
    growth_value = simulate_minimum_essential(model, medium_dict, minimum, anaerobic)
    doubling_time = (np.log(2)/growth_value)*60 if growth_value != 0 else 0
    medium_name = medium['medium'][0] if not anaerobic else f'{medium["medium"][0]}[-O2]'
    exchanges = [[medium_name], minimum,
                 missing_exchanges, [growth_value], [doubling_time]]
    df_growth = pd.DataFrame(exchanges,
                             ['medium',
                              'essential',
                              'missing exchanges',
                              'growth_value',
                              'doubling_time [min]']).T
    return df_growth

# recheck
def get_growth_selected_media(model: cobraModel, media: list[str], basis: str, anaerobic: bool) -> pd.DataFrame:
    """Simulates growth on all given media

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - media (list[str]): Ids of media to simulate on
        - basis (str): Either default_uptake (adding metabs from default) or minimal_uptake (adding metabs from minimal medium)
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        pd.DataFrame: Information on growth behaviour on given media
    """
    growth = pd.DataFrame()
    for medium in media:
        medium_df = load_medium_from_db_for_growth(medium)
        if (basis == 'default_uptake'):
            growth_one = growth_one_medium_from_default(model, medium_df, anaerobic)
        elif (basis == 'minimal_uptake'):
            growth_one = growth_one_medium_from_minimal(model, medium_df, anaerobic)
        growth = pd.concat([growth, growth_one], ignore_index=True)
    return growth

# single knockout simulation
# NOT SIMPLY ESSENTIALS
# entkoppelt
def get_essential_reactions(model: cobraModel) -> list[str]:
    """Knocks out each reaction, if no growth is detected the reaction is seen as essential

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: BiGG Ids of essential reactions
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

    return ess

# same issue as above
# entkoppelt
def get_essential_reactions_via_bounds(model: cobraModel) -> list[str]:
    """Knocks out reactions by setting their bounds to 0, if no growth is detected the reaction is seen as essential


    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: BiGG Ids of essential reactions
    """
    medium = model.medium
    ess = []
    for content in medium.keys():
        model.reactions.get_by_id(content).lower_bound = 0.0
        solution = model.optimize().objective_value
        if solution < 1e-9:
            ess.append(content)
        model.reactions.get_by_id(content).lower_bound = -10.0

    return ess

# rename: find_growth_enhancing_exchanges
def find_additives(model:cobraModel, base_medium: dict) -> pd.DataFrame:
    """Iterates through all exchanges to find metabolites that lead to a higher growth rate compared to the growth rate yielded on the base_medium

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - base_medium (dict): Exchanges as keys and their flux bound as value (f.ex {'EX_glc__D_e' : 10.0})

    Returns:
        pd.DataFrame: Exchanges sorted from highest to lowest growth rate improvement
    """
    with model:
        medium = model.medium
        model.medium = base_medium
        sol = model.optimize()
    base_growth = sol.objective_value
    print(base_growth)

    enhancement = {}
    for ex in list(model.exchanges):
        with model:
            medium = model.medium
            base_medium[ex.id] = 10.0
            model.medium = base_medium
            sol = model.optimize()
            if sol.objective_value > base_growth:
                enhancement[ex.id] = sol.objective_value - base_growth
            base_medium.pop(ex.id)

    adds = pd.DataFrame(enhancement.items(), columns=[
                        'exchange', 'diff']).sort_values(by=['diff'], ascending=False)

    return adds
