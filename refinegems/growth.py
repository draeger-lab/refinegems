#!/usr/bin/env python
""" Provides functions to simulate growth on any medium

Tailored to work with the media denoted in the local db, should work with any medium as long as its defined in a csv with ; as delimiter and BiGG Ids for the compounds. Use refinegems.io.load_medium_custom and hand this to the growth_one_medium_from_default or growth_one_medium_from_minimum function.
"""

import logging
import pandas as pd
import numpy as np
from refinegems.io import load_medium_from_db
from cobra.medium import minimal_medium
from cobra import Reaction
from cobra import Model as cobraModel

__author__ = "Famke Baeuerle"


def set_fluxes_to_simulate(reaction: Reaction) -> Reaction:
    """Helper function: Set flux bounds to -1000.0 and 1000.0 to enable model simulation with growth_one_medium_from_minimal/default

    Args:
        - reaction (Reaction): Reaction with unusable flux bounds

    Returns:
        Reaction: Reaction with usable flux bounds
    """
    reaction.bounds = (-1000.0, 1000.0)
    return reaction


def get_default_uptake(model: cobraModel) -> list[str]:
    """Determines which metabolites are used in the standard medium

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: Metabolites consumed in standard medium
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


def get_minimal_uptake(model: cobraModel) -> list[str]:
    """Determines which metabolites are used in a minimal medium

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: Metabolites consumed in minimal medium
    """
    with model:
        minimal = minimal_medium(model)
    return list(minimal.index)


def get_default_secretion(model: cobraModel) -> list[str]:
    """Checks fluxes after FBA, if positive the metabolite is produced

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: BiGG Ids of produced metabolites
    """
    with model:
        fluxes = model.optimize().fluxes
        default_secretion = []
        for index, value in fluxes.items():
            if value > 0:
                default_secretion.append(index)
    return default_secretion


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
            for reaction in model.reactions:
                set_fluxes_to_simulate(reaction)
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
            for reaction in model.reactions:
                set_fluxes_to_simulate(reaction)
            model.medium = new_medium
        sol = model.optimize()
    return sol.objective_value


def get_all_minimum_essential(model: cobraModel, media: list[str]) -> pd.DataFrame:
    """Returns metabolites necessary for growth and not in media

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - media (list[str]): Containing the names of all media for which the growth essential metabolites not contained in the media should be returned

    Returns:
        pd.DataFrame: information on different media which metabs are missing
    """
    default_uptake = get_default_uptake(model)
    mins = pd.DataFrame()
    for medium in media:
        medium_df = load_medium_from_db(medium)
        missing_exchanges = get_missing_exchanges(model, medium_df)
        medium_dict = modify_medium(medium_df, missing_exchanges)
        essential = find_missing_essential(model, medium_dict, default_uptake)
        minimum = find_minimum_essential(medium_df, essential)
        mins[medium['medium'][0]] = pd.Series(minimum)
    return mins


def growth_one_medium_from_default(model: cobraModel, medium: pd.DataFrame, anaerobic: bool) -> pd.DataFrame:
    """Simulates growth on given medium, adding missing metabolites from the default uptake

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - medium (pd.DataFrame): Dataframe with medium definition
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        pd.DataFrame: Information on growth behaviour on given medium
    """
    default_uptake = get_default_uptake(model)
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
    minimal_uptake = get_minimal_uptake(
        model)  # use this instead of default_uptake
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
        medium_df = load_medium_from_db(medium)
        if (basis == 'default_uptake'):
            growth_one = growth_one_medium_from_default(model, medium_df, anaerobic)
        elif (basis == 'minimal_uptake'):
            growth_one = growth_one_medium_from_minimal(model, medium_df, anaerobic)
        growth = pd.concat([growth, growth_one], ignore_index=True)
    return growth


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
