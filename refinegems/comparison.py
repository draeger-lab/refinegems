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
from libsbml import Model as libModel
from cobra import Model as cobraModel
from refinegems.io import load_medium_from_db, search_sbo_label
from refinegems.growth import growth_one_medium_from_default, growth_one_medium_from_minimal
from refinegems.investigate import initial_analysis, get_reactions_per_sbo

__author__ = "Famke Baeuerle"

def plot_initial_analysis(models: list[libModel]):
    """Creates bar plot of number of entities per Model

    Args:
        models (list[libModel]): Models loaded with libSBML

    Returns:
        plot: pandas plot object
    """
    numbers = pd.DataFrame([initial_analysis(model) for model in models], columns=['model', 'metabolites', 'reactions', 'genes'])
    ax = numbers.set_index('model').plot.bar(y=['metabolites', 'reactions', 'genes'], figsize=(8, 5), cmap='Paired', rot=0)
    # commented is possibility to integrate memote scores
    #numbers.set_index('model').plot(y='Memote score', ax=ax, use_index=False, linestyle=':', secondary_y='Memote score', color='k', marker='D', legend=True)
    #ax.right_ax.set_ylabel('Memote score [%]')
    #ax.right_ax.legend(loc='upper right', bbox_to_anchor=[0.98, 0.9])
    #ax.right_ax.set_ylim([75, 95])
    ax.legend(title=False, loc='upper left', ncol=3, frameon=False)
    ylim = numbers.drop('model', axis=1).max().max() + 200
    ax.set_ylim([0,ylim])
    ax.set_xlabel('')
    ax.tick_params(axis='x',which='both', bottom=False,top=False)
    return ax

def get_sbo_mapping_multiple(models: list[libModel]) -> pd.DataFrame:
    """Determines number of reactions per SBO Term and adss label of SBO Terms

    Args:
        models (list[libModel]): Models loaded with libSBML

    Returns:
        pd.DataFrame: SBO Terms, no of reactions per model and SBO Label
    """
    mappings = {}
    for model in models:
        mappings[model.id] = get_reactions_per_sbo(model)
    df = pd.DataFrame.from_dict(mappings)
    df = df.reset_index().rename({'index': 'SBO-Term'}, axis=1)
    df['SBO-Name'] = df['SBO-Term'].apply(search_sbo_label)
    return df

def plot_rea_sbo_multiple(models: list[libModel], rename=None):
    """Plots reactions per SBO Term in horizontal bar chart with stacked bars for the models

    Args:
        models (list[libModel]): Models loaded with libSBML
        rename (dict, optional): Rename model ids to custom names. Defaults to None.

    Returns:
        plot: pandas plot object
    """
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

def plot_venn(models: list[cobraModel], entity: str, perc: bool=False):
    """Creates venn diagram to show the overlap of model entities

    Args:
        models (list[cobraModel]): Models loaded with cobrapy
        entity (str): Compare on metabolite|reaction
        perc (bool, optional): True if percentages should be used. Defaults to False.

    Returns:
        plot: venn diagram
    """
    intersec = {}
    for model in models:
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

def plot_heatmap_dt(growth: pd.DataFrame):
    """Creates heatmap of simulated doubling times with additives
    
    Args:
        growth (pd.DataFrame): Containing growth data from simulate_all
        
    Returns:
        plot: sns heatmap plot
    """
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

def plot_heatmap_binary(growth: pd.DataFrame):
    """Creates a plot were if growth without additives is possible is marked blue otherwise red

    Args:
        growth (pd.DataFrame): Containing growth data from simulate_all
        
    Returns:
        plot: sns heatmap plot
    """
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

def simulate_all(models: list[cobraModel], media: list[str], basis: str) -> pd.DataFrame:
    """Does a run of growth simulation for multiple models on different media

    Args:
        models (list[cobraModel]): Models loaded with cobrapy
        mediumpath (string): Path to csv containing medium definitions
        media (list): Media of interest (f.ex. LB, M9, ...)
        basis (string): Either default_uptake (adding metabs from default) or minimal_uptake (adding metabs from minimal medium)

    Returns:
        pd.DataFrame: table containing the results of the growth simulation
    """
    growth = pd.DataFrame()
    for medium_id in tqdm(media):
        medium = load_medium_from_db(medium_id)
        for model in models:
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

