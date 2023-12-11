#!/usr/bin/env python
""" Provides functions to simulate growth on any medium and other functionalities replated to growth.
"""

import logging
import pandas as pd
import numpy as np
# from refinegems.io import load_medium_from_db_for_growth # only needed for the old ones
from refinegems.database import medium
from refinegems.io import load_multiple_models, load_model_cobra
from refinegems import reports
from cobra import Model as cobraModel
import cobra
import re
from typing import Literal
import yaml
import matplotlib.pyplot as plt
import warnings

__author__ = "Famke Baeuerle and Carolin Brune"

############################################################################
# functions
############################################################################

# @TEST
def set_bounds_to_default(model: cobraModel, reac_bounds:None|str|tuple[float] = None):
    """Set the reactions bounds of a model to given default values.
    (Ir)reversibility is retained.

    Args:
        model (cobraModel): The model loaded with COBRApy.
        reac_bounds (None|str|tuple[float], optional): The setting for the new reaction bounds. 
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
def get_uptake(model: cobraModel, type: str, exchange_regex:str='^EX') -> list[str]:
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
# @WARNING 
def find_growth_essential_exchanges(model: cobraModel, growth_medium: dict, standard_uptake: list[str]|None) -> list[str]:
    """Find exchanges in a medium (with or without supplements) essential for the growth.
    @WARNING only tests single deletions currently

    Args:
        model (cobraModel): The model to be tested.
        growth_medium (dict): The medium in dictionary form, ready to be added to the model.
        standard_uptake (list[str]|None): Option to add a second medium list as supplements.

    Returns:
        list[str]: The list of exchanges essential for growth.
    """
    with model:
        if not standard_uptake:
            # combine standard and second medium, set all fluxes to 10.0 for standard
            standard_medium = {i: 10.0 for i in standard_uptake}
            new_medium = {**growth_medium, **standard_medium}
        else:
            new_medium = growth_medium
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

    return essential


# @TEST
def find_additives_to_enable_growth(model: cobraModel, growth_medium: dict, standard_uptake: list[str], combine:bool=False):
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
    
    # find essential exchange reactions
    essential = find_growth_essential_exchanges(model, growth_medium, standard_uptake)

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
# @RENAMED 
# @RESTRUCTURED
def get_metabs_essential_for_growth_wrapper(model: cobraModel, media: list[medium.Medium], only_additives:bool=True) -> dict:
    """
    Returns metabolites necessary for growth and not in media

    Args:
        - model (cobraModel): Model loaded with COBRApy
        - media (list[medium.Medium]): Containing all media for which the growth essential metabolites not contained in the media should be returned
        - only_additives(bool, optional): Flag to only return the supplemented exchanges (True) or all essential ones (False).
            Defaults to True.

    Returns:
        dict: information on different media which metabs are missing (key: name of medium, values: list of exchanges)
    """
    
    default_uptake = get_uptake(model,'std')
    ess = {}
    for medium in media:
        
        # convert to a cobrapy medium
        exported_m = medium.medium_to_model(medium=medium, model=model, add=False)
        # get essentials
        essential =  find_growth_essential_exchanges(model, exported_m, default_uptake)
        # check for only_additives flag
        if only_additives:
            essential = [_ for _ in essential if _ not in exported_m.keys()]
        # add to dict
        ess[medium.name] = essential

    return ess


# @TEST
def growth_sim_single(model: cobraModel, m: medium.Medium, supplement:Literal[None,'std','min'] = None) -> reports.SingleGrowthSimulationReport:
    """Simulate the growth of a model on a given medium.

    Args:
        model (cobraModel): The model.
        m (medium.Medium): The medium.
        supplement (Literal[None,'std','min'], optional): Flag to add additvites to the model to ensure growth. Defaults to None (no supplements).
            Further options include 'std' for standard uptake and 'min' for minimal uptake supplementation.

    Returns:
        reports.SingleGrowthSimulationReport: Object with the simulation results
    """

    with model:

        # convert to a cobrapy medium
        exported_m = medium.medium_to_model(medium=m, model=model, 
                                            namespace='BiGG', 
                                            default_flux=10.0, 
                                            replace=False, 
                                            double_o2=False, 
                                            add=False)

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
        report = reports.SingleGrowthSimulationReport(model_name = model.id, medium_name = m.name)
        report.growth_value = model.optimize().objective_value
        report.doubling_time = (np.log(2)/report.growth_value)*60 if report.growth_value != 0 else 0
        report.additives = [_ for _ in new_m if _ not in exported_m]
        report.no_exchange = [_ for _ in m.export_to_cobra().keys() if _ not in exported_m]

    return report


def growth_sim_multi(models: cobraModel|list[cobraModel], media: medium.Medium|list[medium.Medium], supplement_modes:list[Literal['None','min','std']]|None|Literal['None','min','std']=None) -> reports.GrowthSimulationReport:
    """Simulate the growth of (at least one) models on (at least one) media.

    Args:
        models (cobraModel | list[cobraModel]): A COBRApy model or a list of multiple.
        media (medium.Medium | list[medium.Medium]): A refinegems Medium object or a list of multiple.
        supplement_modes (list[Literal[None,'min','std']] | None | Literal[None, 'min', 'std'], optional): Option to supplement the media to enable growth.
            Default to None. Further options include a list with one entry for each medium or a string to set the same default for all.
            The string can be 'min', 'std' or None.

    Returns:
        reports.GrowthSimulationReport: The compiled information of the simulation results.
    """

    # check input 
    if type(models) != list:
        models = [models]
    if type(media) != list:
        media = [media]
    if type(supplement_modes) != list:
        supplement_modes = [supplement_modes * len(media)]

    # simulate the growth of the models on the different media
    report = reports.GrowthSimulationReport()
    for mod in models:
        for med,supp in zip(media, supplement_modes):
            r = growth_sim_single(mod, med, supplement=supp)
            report.add_sim_results(r)

    return report     


# @IDEA: more options for fluxes
# @TODO/@IDEA: validity check nefore parsing
def read_media_config(yaml_path:str):

    media_list = []
    supplement_list = []

    with open(yaml_path, 'r') as stream:

        loaded = yaml.safe_load(stream)

        # ........................
        # @TODO / MAYBE
        # check validity of input?
        # ........................

        params = loaded['params'] if 'params' in loaded.keys() else None
        media = loaded['media'] if 'media' in loaded.keys() else None

        # handle internal/in-build media compositions
        if media:

            for name,p in media.items():

                # load base medium from either  
                # base
                if p and 'base' in p.keys():
                    new_medium = medium.load_medium_from_db(p['base'])
                    new_medium.name = name
                # extern
                elif p and 'external_base' in p.keys():
                    new_medium = medium.read_external_medium(p['external_base'])
                    new_medium.name = name
                # default name
                else:
                    new_medium = medium.load_medium_from_db(name)
                
                # add additional media 
                # from in-build database
                if p and 'add' in p.keys():
                    for a in p['add']:
                        if a in medium.SUBSET_MEDIA_MAPPING.keys():
                            new_medium = new_medium.add_subset(a)
                        else:
                            new_medium = new_medium + medium.load_medium_from_db(a)
                # from external
                if p and 'add_external' in p.keys():
                    for a in p['add_external']:
                        new_medium = new_medium + medium.read_external_medium(a)
                
                # check anaerobic / aerobic settings
                if p and 'aerobic' in p.keys():
                    if p['aerobic']:
                        new_medium.make_aerobic()
                    else:
                        new_medium.make_anaerobic()
                else:
                    if params and 'aerobic' in params.keys():
                        if params['aerobic']:
                            new_medium.make_aerobic()
                        else:
                            new_medium.make_anaerobic()

                # set default flux
                if p and 'default_flux' in p.keys():
                    new_medium.set_default_flux(p['default_flux'], replace=True)
                elif params and 'default_flux' in params.keys():
                    new_medium.set_default_flux(params['default_flux'], replace=True)

                # set o2_percentage
                if p and 'o2_percent' in p.keys():
                    new_medium.set_oxygen_percentage(p['o2_percent'])
                elif params and 'o2_percemt' in params.keys():
                    new_medium.set_oxygen_percentage(params['o2_percent'])

                # supplement settings
                if p and 'supplement' in p.keys():
                    supplement_list.append(p['supplement'])
                elif params and 'supplement' in params:
                    supplement_list.append(params['supplement'])
                else:
                    supplement_list.append(None)
                    
                # append medium to list
                media_list.append(new_medium)

    return (media_list,supplement_list)


# @IDEA : choose different namespaces for the media
def growth_analysis(models:cobra.Model|str|list[str]|list[cobra.Model],
                    media:medium.Medium|list[medium.Medium]|str,
                    supplements:None|list[Literal[None,'std','min']]=None,
                    retrieve:Literal['report','plot','both']='plot') -> reports.GrowthSimulationReport|plt.Figure|tuple:

    # read-in all models into list
    # ----------------------------
    mod_list = []
    match models:
        case list():
            # list as input
            if len(models) > 0:
                # if list entries are paths
                if all(isinstance(_, str) for _ in models):
                    mod_list = load_multiple_models(models, package='cobra')
                # if list entries are already cobra.Models
                elif all(isinstance(_, cobra.Model) for _ in models):
                    mod_list = models
                # @TODO
                # option for mixed list?
                else:
                    raise TypeError('Unknown or mixed types in model list.')
            else:
                raise KeyError('Empty list for models detected.')
        # single model as input
        case cobra.Model():
            mod_list = [cobra.Model]
        # single string as input
        case str():
            mod_list = [load_model_cobra(models)]
        # unknown input
        case _:
            raise ValueError(F'Unknown input type for models: {type(models)}')
        
    # collect all media into list
    # ---------------------------
    media_list = []
    match media:
        # single medium
        case medium.Medium():
            media_list = [media]
        # list of media
        case list():
            if all(isinstance(_,medium.Medium) for _ in media):
                media_list = media
            else:
                raise TypeError('Unknown type found in media, should be list fo medium.Medium.')
        # string - connection to YAML config file
        case str():
            media, supplements = read_media_config(media)
        # unknown input
        case _:
            raise ValueError(f'Unknown input for media: {media}')

    # run simulation
    # --------------
    report = growth_sim_multi(mod_list, media, supplements)

    # save / visualise report 
    # -----------------------
    match retrieve:
        case 'report':
            return report
        case 'plot':
            return report.plot_growth()
        case 'both':
            return (report, report.plot_growth)
        case _:
            raise ValueError(f'Unknown input for retrieve: {retrieve}')


# @RENAMED
def get_essential_reactions_via_single_knockout(model: cobraModel) -> list[str]:
    """Knocks out each reaction, if no growth is detected the reaction is seen as essential

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: Ids of essential reactions
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

# @RENAMED
def get_essential_exchanges_via_bounds(model: cobraModel) -> list[str]:
    """Knocks out reactions by setting their bounds to 0, if no growth is detected the reaction is seen as essential


    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        list[str]: Ids of essential reactions
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

# @RENAMED
def find_growth_enhancing_exchanges(model:cobraModel, base_medium: dict) -> pd.DataFrame:
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


# @TEST
# from SPECIMEN
# auxothrophy test 
def test_auxotrophies(model:cobraModel, media_list:list[medium.Medium], namespace:Literal['BiGG']='BiGG') -> dict[dict[str,float]]:
    """Test for amino acid auxothrophies for a model and a list of media.

    Tests, if the model growths on the media and if and with what fluxes the
    20 proteinogenic amino acids are produced by temporarily adding a
    sink reaction for each of the amino acids to the model as the objective function.

    Args:
        model (cobraModel): The model to be tested. Loaded with COBRApy.
        media_list (list[medium.Medium]): List of media to be tested.
        namespace (Literal['BiGG'], optional): String for the namespace to be used for the model. 
            Current options include 'BiGG'.
            Defaults to 'BiGG'.

    Raises:
        ValueError: Unknown input for namespace parameter.

    Returns:
        dict[dict[str,float]]: The test results as a dictionary of the media names as keys with dictionaries of the amino acids names and growth flux rates of the sink reaction as values.
    """

    results = {}

    # get amino acids from database
    amino_acids = medium.Medium('aa').add_subset(type='aa')
    aa_list = set(amino_acids.substance_table['name'])

    # iterate over all media
    for med in media_list:
        auxotrophies = {}
        # then iterate over all amino acids
        for a in aa_list:

            entry = med.substance_table[(med.substance_table['name'] == a) & (med.substance_table['db_type'] == namespace)]

            with model as m:

                # first set the medium
                medium.medium_to_model(m, med, namespace=namespace, default_flux=10.0, replace=False, double_o2=False, add=True)
                
                # check namespace availability
                if len(entry) == 0:
                    warnings.warn('Amino acid {a} has no identifier for your chosen namespace {namespace}. Please contact support if you want to add one.')
                    auxotrophies[a] = m.optimize().objective_value
                else:
                    
                    # create and check IDs for the chosen namespace
                    match namespace:
                        case 'BiGG':
                            internal_meta = ""
                            for np_id in entry['db_id']:
                                if np_id + '_c' in [_.id for _ in m.metabolites]:
                                    internal_meta = np_id + '_c'
                                    break
                            if internal_meta == "":
                                warnings.warn(F'No identifier matched in cytosol for {a}.')

                            exchange_reac = 'EX_' + internal_meta
                            sink_reac = F'sink_{internal_meta}_tmp'
                        case _:
                            raise ValueError('Unknown namespace: {namespace}. Cannot create IDs.')
                                
                # create a pseudo reaction -> a sink reaction for the amino acid
                # to use as the new objective
                if internal_meta == "":
                    pass
                else:
                    m.add_boundary(m.metabolites.get_by_id(internal_meta), type='sink', reaction_id=sink_reac)
                    m.objective = sink_reac
                    # if existent, close the exchange reaction
                    if exchange_reac in [_.id for _ in m.exchanges]:
                        m.reactions.get_by_id(exchange_reac).lower_bound = 0.0
                        m.reactions.get_by_id(exchange_reac).upper_bound = 0.0
                # and calculate the new objective
                auxotrophies[a] = m.optimize().objective_value()

            # add the current test results to the list of all results
            results[med.name] = auxotrophies

    return results

