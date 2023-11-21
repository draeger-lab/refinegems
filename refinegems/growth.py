#!/usr/bin/env python
""" Provides functions to simulate growth on any medium

Tailored to work with the media denoted in the local db, should work with any medium as long as its defined in a csv with ; as delimiter and BiGG Ids for the compounds. Use refinegems.io.load_medium_custom and hand this to the growth_one_medium_from_default or growth_one_medium_from_minimum function.
"""

import logging
import pandas as pd
import numpy as np
from refinegems.io import load_medium_from_db_for_growth
from refinegems.database import medium
from refinegems import reports
from cobra import Reaction
from cobra import Model as cobraModel
import cobra
import re

__author__ = "Famke Baeuerle and Carolin Brune"

############################################################################
# functions
############################################################################

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


# @TEST
def find_additives_to_enable_growth(model: cobraModel, growth_medium: dict, standard_uptake: list[str], combine=False):
    """Based on a new medium for growth and a standard one the model already growths on, find additives from the standard, 
    which can be added to the new one to enable growths.

    @WARNING this functions currently only check via single deletions not via correlated ones. MIght lead to problems in 
    complecated cases, e.g. missing carbon source but multiple ones in the standard medium.

    Args:
        model (cobraModel): The model to add the medium to.
        growth_medium (dict): A medium definition, ready to be added to a COBRApy model.
            This is for the new medium for growth testing.
        standard_uptake (list[str]): A list of exchange reactions, e.g. Output of get_uptake.
            This is for the old medium where the model can grow on.
        combine (bool, optional): Flag to directly combine the additives with the new medium or just return the additive's IDs.
            Defaults to False (returns reaction IDs).

    Returns:
        list[str] or dict: List of the exchange reaction IDs of the additives or the supplemented medium, if combine is set to True.
    """
    
    with model:
        # combine standard and second medium, set all fluxes to 10.0 for standard
        standard_medium = {i: 10.0 for i in standard_uptake}
        new_medium = {**growth_medium, **standard_medium}
        # add new medium to model
        try:
            model.medium = new_medium
        except(ValueError):
            logging.info('Change upper bounds to COBRApy defaults to make model simulatable.')
            set_bounds_to_default(model)
            model.medium = new_medium
        # find essentials 
        # @WARNING: only single deletions are tested
        essential = []
        for metab in new_medium.keys():
            with model:
                model.reactions.get_by_id(metab).lower_bound = 0
                sol = model.optimize()
                if sol.objective_value < 1e-5:  # and sol.objective_value > -1e-9: # == 0 no negative growth!
                    essential.append(metab)

    # find the essential compounds not in the growth medium
    additives = []
    for metab in essential:
        if metab not in growth_medium.keys():
            additives.append(metab)

    # return ... 
    if combine:
        # ... the supplemented medium
        supplements = {_: 10.0 for _ in additives}
        suppl_medium = {**growth_medium,**supplements}
        return suppl_medium
    else:
        # ... the list of suplements
        return additives

# @TEST
def growth_sim_single(model: cobraModel, m: medium.Medium, supplement = None, anaerobic = False) -> reports.SingleGrowthSimulationReport:
    """Simulate the growth of a model on a given medium.

    Args:
        model (cobraModel): The model.
        m (medium.Medium): The medium.
        supplement (str, optional): Flag to add additvites to the model to ensure growth. Defaults to None (no supplements).
            Further options include 'std' for standard uptake and 'min' for minimal uptake supplementation.
        anaerobic (bool, optional): Set medium to anaerobic/aerobic. Defaults to False (aerobic growth).

    Returns:
        reports.SingleGrowthSimulationReport: Object with the simulation results
    """

    with model:
        # check for (an)aerobic conditions
        if anaerobic and m.is_aerobic():
            m.make_anaerobic()
        if not anaerobic and not m.is_aerobic():
            m.make_aerobic(flux=10.0)

        # convert to a cobrapy medium
        exported_m = medium.medium_to_model(medium=m, model=model, add=False)

        # supplement, if tag is set
        match supplement:
            case 'std':
                uptake = get_uptake(model,'std')
                new_m = find_additives_to_enable_growth(model, exported_m, uptake, combine=True)
            case 'min':
                uptake = get_uptake(model,'min')
                new_m = find_additives_to_enable_growth(model, exported_m, uptake, combine=True)
            case _:
                new_m = exported_m
    
        # add medium to model 
        try:
            model.medium = new_m
        except(ValueError):
            logging.info('Change upper bounds to 1000.0 and lower bounds to -1000.0 to make model simulatable.')
            set_bounds_to_default(model)
            model.medium = new_m

        # simulate growth
        report = reports.SingleGrowthSimulationReport(model_name = model.name, medium_name = m.name)
        report.growth_value = model.optimize().objective_value
        report.doubling_time = (np.log(2)/report.growth_value)*60 if report.growth_value != 0 else 0
        report.additives = [_ for _ in new_m if _ not in exported_m]
        report.no_exchange = [_ for _ in m.export_to_cobra().keys() if _ not in exported_m]

    return report


# @TEST
def growth_sim_multi(models: cobraModel|list[cobraModel], media: medium.Medium|list[medium.Medium]) -> reports.GrowthSimulationReport:
    """Simulate the growth of (at least one) models on (at least one) media.

    Args:
        models (cobraModel | list[cobraModel]): A COBRApy model or a list of multiple.
        media (medium.Medium | list[medium.Medium]): A refinegems Medium object or a list of multiple.

    Returns:
        reports.GrowthSimulationReport: The compiled information of the simulation results.
    """

    # check input 
    if type(models) != list:
        models = [models]
    if type(media) != list:
        media = [media]

    # simulate the growth of the models on the different media
    report = reports.GrowthSimulationReport()
    for mod in models:
        for med in media:
            o2_check = med.is_aerobic()
            r = growth_sim_single(mod, med, anaerobic = not o2_check)
            report.add_sim_results(r)

    return report     


# @TODO
# main objective: read in models and media from input (command line, YAML etc.) 
# -> compile a complete list media 
# -> run simulation 
# -> visulise also here?
def growth_analysis():
    pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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


def simulate_all(models: list[cobraModel], media: list[str], basis: str, anaerobic: bool) -> pd.DataFrame:
    """Does a run of growth simulation for multiple models on different media

    Args:
        - models (list[cobraModel]): Models loaded with cobrapy
        - media (list[str]): Media of interest (f.ex. LB, M9, ...)
        - basis (str): Either default_uptake (adding metabs from default) or minimal_uptake (adding metabs from minimal medium)
        - anaerobic (bool): If True 'EX_o2_e' is set to 0.0 to simulate anaerobic conditions

    Returns:
        pd.DataFrame: table containing the results of the growth simulation
    """
    growth = pd.DataFrame()
    for medium_id in tqdm(media):
        medium = load_medium_from_db(medium_id)
        for model in models:
            essentials_given = False
            if (basis=='default_uptake'):
                growth_one = growth_one_medium_from_default(model, medium, anaerobic).drop('missing exchanges', axis=1)
            elif (basis == 'minimal_uptake'):
                growth_one = growth_one_medium_from_minimal(model, medium, anaerobic).drop('missing exchanges', axis=1)
            if growth_one['essential'].dropna().size == 0:
                essentials_given = True
            else:
                growth_list = growth_one['essential'].dropna().to_list()
                growth_string = ', '.join(growth_list)
                essentials_given = growth_string
            growth_one = growth_one.drop('essential', axis=1)
            growth_one['complete'] = essentials_given
            growth_one = growth_one.dropna()
            growth_one['model'] = model.id
            growth_one = growth_one[['model', 'medium', 'doubling_time [min]', 'growth_value', 'complete']]
            growth_one['doubling_time [min]'].astype(float).round(2)
            growth_one['growth_value'].astype(float).round(2)
            growth = growth.append(
                growth_one, 
                ignore_index=True)

    return growth




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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
