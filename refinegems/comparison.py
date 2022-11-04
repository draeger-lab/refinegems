#!/usr/bin/env python
""" Provides functions to compare multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""

import pandas as pd
from tqdm import tqdm
from refinegems.load import load_model_cobra, load_all_media_from_db
from refinegems.growth import growth_one_medium_from_default, growth_one_medium_from_minimal

__author__ = "Famke Baeuerle"


def simulate_all(model_list, mediumpath, media, basis):
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
    selected_media = [x for x in all_media if x['medium'][0] in media]
    for medium in tqdm(selected_media):
        for model_path in model_list:
            model = load_model_cobra(model_path)
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
