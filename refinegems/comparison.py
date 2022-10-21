#!/usr/bin/env python
""" Provides functions to compare multiple models

Can mainly be used to compare growth behaviour of multiple models. All other stats are shown in the memote report.
"""

__author__ = "Famke Baeuerle"
import pandas as pd
from tqdm import tqdm
from refinegems.load import load_model_cobra, load_all_media_from_db
from refinegems.growth import get_growth_one_medium
from refinegems.investigate import get_egc # implement here for multiple models

def simulate_all(model_list, mediumpath, media):
    growth = pd.DataFrame()
    all_media = load_all_media_from_db(mediumpath)
    selected_media = [x for x in all_media if x['medium'][0] in media]
    for medium in tqdm(selected_media):
        for model_path in model_list:
            model = load_model_cobra(model_path)
            essentials_given = False
            growth_one = get_growth_one_medium(model, medium).drop('missing', axis=1)
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
    #print(growth.to_latex(index = False))
    return growth
