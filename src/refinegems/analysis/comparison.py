#!/usr/bin/env python
""" Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""


__author__ = "Famke Baeuerle"

################################################################################
# requirements
################################################################################

import pandas as pd

from cobra import Model as cobraModel
from libsbml import Model as libModel
from venn import venn

from ..utility.io import search_sbo_label
from .investigate import get_reactions_per_sbo

################################################################################
# functions
################################################################################

def get_sbo_mapping_multiple(models: list[libModel]) -> pd.DataFrame:
    """Determines number of reactions per SBO Term and adds label of SBO Terms

    Args:
        - models (list[libModel]): 
            Models loaded with libSBML

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
        - models (list[libModel]): 
            Models loaded with libSBML
        - rename (dict, optional): 
            Rename model ids to custom names. 
            Defaults to None.

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
        - models (list[cobraModel]): 
            Models loaded with cobrapy
        - entity (str): 
            Compare on metabolite|reaction
        - perc (bool, optional): 
            True if percentages should be used. 
            Defaults to False.
        - rename (dict, optional): 
            Rename model ids to custom names. 
            Defaults to None.

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

