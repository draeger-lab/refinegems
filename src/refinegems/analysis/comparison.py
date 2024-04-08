#!/usr/bin/env python
""" Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""


__author__ = "Famke Baeuerle"

################################################################################
# requirements
################################################################################

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm
from venn import venn
from libsbml import Model as libModel
from cobra import Model as cobraModel

from ..utility.io import search_sbo_label
from .investigate import get_reactions_per_sbo

################################################################################
# functions
################################################################################

def get_sbo_mapping_multiple(models: list[libModel]) -> pd.DataFrame:
    """Determines number of reactions per SBO Term and adds label of SBO Terms

    Args:
        - models (list[libModel]): Models loaded with libSBML

    Returns:
        pd.DataFrame: SBO Terms, number of reactions per Model and SBO Label
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
        - models (list[libModel]): Models loaded with libSBML
        - rename (dict, optional): Rename model ids to custom names. Defaults to None.

    Returns:
        plot: Pandas stacked barchart
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

def plot_venn(models: list[cobraModel], entity: str, perc: bool=False, rename=None):
    """Creates Venn diagram to show the overlap of model entities

    Args:
        - models (list[cobraModel]): Models loaded with cobrapy
        - entity (str): Compare on metabolite|reaction
        - perc (bool, optional): True if percentages should be used. Defaults to False.
        - rename (dict, optional): Rename model ids to custom names. Defaults to None.

    Returns:
        plot: Venn diagram 
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
        if rename is not None:
            intersec[rename[model.id]] = set(reas)
        else:
            intersec[model.id] = set(reas)
    if perc:
        fig = venn(intersec, fmt="{percentage:.1f}%")
    else:
        fig = venn(intersec)
    return fig

def plot_heatmap_dt(growth: pd.DataFrame):
    """Creates heatmap of simulated doubling times with additives
    
    Args:
        - growth (pd.DataFrame): Containing growth data from simulate_all
        
    Returns:
        plot: Seaborn Heatmap
    """
    growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
    growth.columns.name=None
    growth.index.names = (None,None)
    growth.index.name=None
    growth.index = growth.index.get_level_values(1)
    growth[growth > 500] = 0
    growth[growth < 0] = 0
    growth.replace([np.inf, -np.inf], 0, inplace=True)
    over_growth = growth.max().max() + 6
    growth.replace(np.nan, over_growth, inplace=True)
    under_growth = growth.min().min() - 5
    vmin= under_growth if under_growth > 1e-5 else 1e-5 #Use same threshhold as in find_missing_essential in growth
    vmax=over_growth - 1
    annot = growth.copy()
    annot = annot.round().astype(int)
    annot[annot < 1e-5] = ''
    annot.replace(over_growth.round().astype(int), 'No data', inplace=True)
    cmap=matplotlib.cm.get_cmap('YlGn').copy()
    cmap.set_under('black')
    cmap.set_over('white')
    fig, ax = plt.subplots(figsize=(10,8))
    sns.heatmap(growth.T, 
                annot=annot.T, 
                annot_kws={"fontsize":15},
                vmin=vmin, 
                vmax=vmax,
                cmap=cmap, 
                linewidth=.5, 
                cbar_kws = {'orientation':'vertical', 'label':'doubling time [min]', 'extend': 'min', 'extendrect':True},
                ax=ax,
                fmt=''
                )
    rotation = 40 if len(growth.index) > 3 else 0
    plt.tick_params(rotation=0, bottom=False, top=False, left=False, right=False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha="right")
    return fig

def plot_heatmap_native(growth: pd.DataFrame):
    """Creates a plot were if growth without additives is possible is marked from yellow to green otherwise black

    Args:
        - growth (pd.DataFrame): Containing growth data from simulate_all
        
    Returns:
        plot: Seaborn Heatmap
    """
    def get_native_growth(row):
        if row['complete'] == True:
            return row['doubling_time [min]']
        else:
            return 0
    
    growth['native_growth'] = growth.apply(get_native_growth, axis=1)
    growth = growth[['medium', 'model', 'native_growth']]
    growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
    growth.columns.name=None
    growth.index.names = (None,None)
    growth.index.name=None
    growth.index = growth.index.get_level_values(1)
    growth[growth > 500] = 0
    growth[growth < 0] = 0
    growth.replace([np.inf, -np.inf], 0, inplace=True)
    over_growth = growth.max().max() + 6
    growth.replace(np.nan, over_growth, inplace=True)
    annot = growth.copy()
    annot = annot.round().astype(int)
    annot[annot == np.nan] = 'No data'
    annot[annot < 1e-5] = ''
    annot.replace(over_growth.round().astype(int), 'No data', inplace=True)
    under_growth = growth.min().min() - 5
    vmin= under_growth if under_growth > 1e-5 else 1e-5 #Use same threshhold as in find_missing_essential in growth
    vmax= over_growth - 1
    cmap=matplotlib.cm.get_cmap('YlGn').copy()
    cmap.set_under('black')
    cmap.set_over('white')
    fig, ax = plt.subplots(figsize=(10,8))
    sns.heatmap(growth.T,
                annot=annot.T, 
                annot_kws={"fontsize":15},
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
                linewidth=.5, 
                ax=ax,
                cbar_kws={'orientation':'vertical', 'label':'doubling time [min]', 'extend': 'min', 'extendrect':True},
                fmt=''
                )
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    rotation = 40 if len(growth.index) > 3 else 0
    plt.tick_params(rotation=0, bottom=False, top=False, left=False, right=False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha="right")
    return fig

