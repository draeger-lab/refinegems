#!/usr/bin/env python
""" Provides functions for adding charges to metabolites

When iterating through all metabolites present in a model, you will find several which have no defined charge (metab.getPlugin('fbc').isSetCharge() = false). This can lead to charge imbalanced reactions. This script takes information on metabolite charges from the ModelSEED database. A charge is automatically added to a metabolite if it has no defined charge and if there is only one charge denoted in ModelSEED. When multiple charges are present, the metabolite and the possible charges are noted and later returned in a dictionary.

It is possible to use the correct_charges_from_db function with other databases. The user just needs to make sure that the compounds dataframe has a 'BiGG' and a 'charge' column.
"""

import pandas as pd
from libsbml import *
from refinegems.load import write_to_file
from refinegems.modelseed import get_modelseed_compounds
import re

__author__ = "Famke Baeuerle"


def correct_charges_from_db(model, compounds):
    """Adds charges taken from given database to metabolites which have no defined charge

    Args:
        model (libsbml-model): model loaded with libsbml
        compounds (df): containing database data with 'BiGG' (BiGG-Ids) and 'charge' (float or int) as columns

    Returns:
        tuple: (model with added charges, metabolites with multiple charges as pd df)
    """
    spe = model.getListOfSpecies()
    mulchar = dict()
    for i in spe:
        if not i.getPlugin('fbc').isSetCharge(
        ):  # we are only interested in metab without charge
            bigg = i.getId()[2:-2]
            if len(compounds[compounds['BiGG']
                   == bigg]['charge'].array) == 1:  # eindeutig
                charge = compounds[compounds['BiGG']
                                             == bigg]['charge'].array[0]
                i.getPlugin('fbc').setCharge(int(charge))
            elif len(compounds[compounds['BiGG'] == bigg]['charge'].array) > 1:
                charges = compounds[compounds['BiGG']
                                              == bigg]['charge'].array
                if all(x == charges[0] for x in charges):
                    charge = charges[0]
                    i.getPlugin('fbc').setCharge(int(charge))
                else:
                    mulchar[bigg] = charges

    return model, mulchar


def correct_charges_modelseed(model, new_file_path, modelseed_path, charge_report_path):
    """wrapper function which completes the steps to charge correction

    Args:
        model (libsbml-model): model loaded with libsbml
        new_file_path (Str): filepath + name for modified model
        modelseed_path (str): path to modelseed compound definition

    Returns:
        dict: BiGG Id and possible charges of metabolites
    """
    modelseed_compounds = get_modelseed_compounds(modelseed_path)
    model_corr, multiple_charges = correct_charges_from_db(
        model, modelseed_compounds)
    write_to_file(model_corr, new_file_path)
    pd.DataFrame.from_dict(multiple_charges, orient='index').to_csv(charge_report_path, sep=',', header=False)
