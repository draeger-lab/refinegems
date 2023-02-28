#!/usr/bin/env python
""" Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm
from venn import venn
from refinegems.io import load_multiple_models, load_medium_from_db, search_sbo_label
from refinegems.growth import growth_one_medium_from_default, growth_one_medium_from_minimal
from refinegems.investigate import initial_analysis, get_reactions_per_sbo

__author__ = "Famke Baeuerle"

def get_sbo_mapping_multiple(models):
    mappings = {}
    for model in models:
        mappings[model.id] = get_reactions_per_sbo(model)
    df = pd.DataFrame.from_dict(mappings)
    df = df.reset_index().rename({'index': 'SBO-Term'}, axis=1)
    df['SBO-Name'] = df['SBO-Term'].apply(search_sbo_label)
    return df

def get_sbo_plot_multiple(model_list: list[str], rename=None):
    models = load_multiple_models(model_list, package='libsbml')
    map = get_sbo_mapping_multiple(models)
    id_list = [mod.id for mod in models]
    map = map[(map[id_list]>3).all(axis=1)]
    map = map.drop('SBO-Term', axis=1).sort_values(id_list[0]).set_index('SBO-Name')
    if rename is not None:
        map = map.rename(rename, axis=1)
    fig = map.plot.barh(stacked=True, width=.8, figsize=(8,10))
    fig.set_ylabel('')
    fig.set_xlabel('number of reactions', fontsize=16)
    fig.legend(loc='lower right')
    return fig

def create_venn(model_list: list[str], entity: str, perc: bool=False): 
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

def create_heatmap(growth: pd.DataFrame):
    growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
    growth.columns.name=None
    growth.index.names = (None,None)
    growth.index.name=None
    growth.index = growth.index.get_level_values(1)
    growth[growth > 500] = np.nan
    growth[growth < 0] = np.nan
    growth.replace([np.inf, -np.inf], np.nan, inplace=True)
    vmin=growth.min().min() - 5
    vmax=growth.max().max() + 5
    fig, ax = plt.subplots(figsize=(10,8))
    sns.heatmap(growth.T, 
                annot=True, 
                annot_kws={"fontsize":15},
                vmin=vmin, 
                vmax=vmax,
                #cmap='crest', 
                linewidth=.5, 
                cbar_kws = {'orientation':'horizontal', 'label':'doubling time [min]'},
                ax=ax,
                fmt='.0f'
                )
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        )
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,
        )
    return fig

def create_binary_heatmap(growth: pd.DataFrame):
    def get_native_growth(row):
        if row == True:
            return 1
        else:
            return 0
    growth['native_growth'] = growth['complete'].apply(get_native_growth)
    growth = growth[['medium', 'model', 'native_growth']]
    growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
    growth.columns.name=None
    growth.index.names = (None,None)
    growth.index.name=None
    growth.index = growth.index.get_level_values(1)
    fig, ax = plt.subplots(figsize=(10,8))
    sns.heatmap(growth.T, 
                annot_kws={"fontsize":15},
                cmap='RdYlBu', 
                linewidth=.5, 
                ax=ax,
                cbar=False,
                )
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        )
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,
        )
    return fig

def simulate_all(model_list: list[str], media: list[str], basis: str) -> pd.DataFrame:
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
    all_models = load_multiple_models(model_list, package='cobra')
    for medium_id in tqdm(media):
        medium = load_medium_from_db(medium_id)
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

