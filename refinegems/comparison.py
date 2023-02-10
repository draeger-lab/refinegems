#!/usr/bin/env python
""" Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""

import pandas as pd
import matplotlib.pyplot as plt
from cobra import Model
from tqdm import tqdm
from venn import venn
from refinegems.load import load_multiple_models, load_all_media_from_db
from refinegems.growth import growth_one_medium_from_default, growth_one_medium_from_minimal
from refinegems.investigate import initial_analysis

__author__ = "Famke Baeuerle"

sbo_mapping={'658': 'passive transport', 
            '176': 'biochemical reaction', 
            '167': 'biochemical or transport reaction',
            '402': 'transfer of a chemical group', 
            '659': 'symporter-mediated transport', 
            '200': 'redox reaction', 
            '233': 'hydroxylation',
            '399': 'decarboxylation', 
            '178': 'cleavage', 
            '403': 'transamination', 
            '215': 'acetylation', 
            '377': 'isomerisation', 
            '657': 'active transport', 
            '216': 'phosphorylation', 
            '401': 'deamination', 
            '376': 'hydrolysis', 
            '217': 'glycosylation', 
            '660': 'antiporter-mediated transport', 
            '654': 'co-transport reaction', 
            '214': 'methylation', 
            '655': 'transport reaction', 
            '627': 'exchange reaction', 
            '632': 'sink reaction', 
            '629': 'biomass production',
            '630': 'ATP maintenance'}

def create_venn(model_list: list[str], entity: str, perc: bool=False) -> plt.figure: 
    all_models = load_multiple_models(model_list, package='cobra')
    intersec = {}
    for model in all_models:
        reas = []
        if entity == 'metabolite':
            for rea in model.metabolites:
                reas.append(rea.id)
        if entity == 'reaction':
            for rea in model.reactions:
                reas.append(rea.id)
        intersec[model.id] = set(reas)
    if perc:
        fig = venn(intersec, fmt="{percentage:.0f}%")
    else:
        fig = venn(intersec)
    return fig

def simulate_all(model_list: list[str], mediumpath: str, media: list[str], basis: str) -> pd.DataFrame:
    """does a run of growth simulation for multiple models on different media

    Args:
        model_list (list): paths to the models of interest (xml files)
        mediumpath (string): path to csv containing medium definitions
        media (list): media of interest (f.ex. LB, M9, ...)
        basis (string): either default_uptake (adding metabs from default) or minimal_uptake (adding metabs from minimal medium)

    Returns:
        df: table containing the results of the growth simulation
    """
    growth = pd.DataFrame()
    all_media = load_all_media_from_db(mediumpath)
    all_models = load_multiple_models(model_list, package='cobra')
    selected_media = [x for x in all_media if x['medium'][0] in media]
    for medium in tqdm(selected_media):
        for model in all_models:
            essentials_given = False
            if (basis=='default_uptake'):
                growth_one = growth_one_medium_from_default(model, medium).drop('missing exchanges', axis=1)
            elif (basis == 'minimal_uptake'):
                growth_one = growth_one_medium_from_minimal(model, medium).drop('missing exchanges', axis=1)
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

