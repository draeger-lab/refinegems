#!/usr/bin/env python
""" Provides functions for adding charges to metabolites

When iterating thorugh all metabolites present in a model, you will find several which have no defined charge (metab.getPlugin('fbc').isSetCharge() = false). This can lead to charge imbalanced reactions. This script takes information on metabolite charges from the ModelSEED database. A charge is automatically added to a metabolite if it has no defined charge and if there is only one charge denoted in ModelSEED. When multiple charges are present, the metabolite and the possible charges are noted and later returned in a dictionary.
"""

import pandas as pd
from libsbml import *
from refinegems.load import write_to_file
import re

__author__ = "Famke Baeuerle"

def get_modelseed_compounds(path):
    """extracts compounds from modelseed which have BiGG Ids

    Args:
        path (str): path to modelseed compound definition

    Returns:
        df: table containing modelseed data
    """
    com = pd.read_csv(path, sep='\t')

    def get_bigg_ids(aliases):
        try:
            aliases_list = aliases.split('|')
            bigg = [x[6:] for x in aliases_list if re.search('BiGG: .*', x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com['BiGG'] = com.apply(lambda row: get_bigg_ids(row['aliases']), axis=1)

    return com.loc[:, ['id', 'name', 'formula', 'mass', 'charge', 'BiGG']].dropna(subset=['BiGG'])

def correct_charges_modelseed(model, modelseed_compounds):
    """Adds charges taken from the ModelSEED database to metabolites which have no defined charge

    Args:
        model (libsbml-model): model loaded with libsbml
        modelseed_compounds (df): containing modelseed data

    Returns:
        tuple: (model with added charges, metabolites with multiple charges)
    """
    spe = model.getListOfSpecies()
    mulchar = dict()
    for i in spe:
        if not i.getPlugin('fbc').isSetCharge(): # we are only interested in metab without charge
            bigg = i.getId()[2:-2]
            if len(modelseed_compounds[modelseed_compounds['BiGG'] == bigg]['charge'].array) == 1: #eindeutig
                charge = modelseed_compounds[modelseed_compounds['BiGG'] == bigg]['charge'].array[0]
                i.getPlugin('fbc').setCharge(int(charge))
            elif len(modelseed_compounds[modelseed_compounds['BiGG'] == bigg]['charge'].array) > 1:
                charges = modelseed_compounds[modelseed_compounds['BiGG'] == bigg]['charge'].array
                if all(x==charges[0] for x in charges):
                    charge = charges[0]
                    i.getPlugin('fbc').setCharge(int(charge))
                else:
                    mulchar[bigg] = charges
    
    return model, mulchar

    
def charges(model, new_file_path, modelseed_path):
    """wrapper function which completes the steps to charge correction 

    Args:
        model (libsbml-model): model loaded with libsbml
        new_file_path (Str): filepath + name for modified model
        modelseed_path (str): path to modelseed compound definition

    Returns:
        dict: BiGG Id and possible charges of metabolites
    """
    modelseed_compounds = get_modelseed_compounds(modelseed_path)
    model_corr, multiple_charges = correct_charges_modelseed(model, modelseed_compounds)
    write_to_file(model_corr, new_file_path)
    
    return multiple_charges