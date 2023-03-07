#!/usr/bin/env python
""" Provides functions for adding charges to metabolites

When iterating through all metabolites present in a model, you will find several which have no defined charge (metab.getPlugin('fbc').isSetCharge() = false). This can lead to charge imbalanced reactions. This script takes information on metabolite charges from the ModelSEED database. A charge is automatically added to a metabolite if it has no defined charge and if there is only one charge denoted in ModelSEED. When multiple charges are present, the metabolite and the possible charges are noted and later returned in a dictionary.

It is possible to use the correct_charges_from_db function with other databases. The user just needs to make sure that the compounds dataframe has a 'BiGG' and a 'charge' column.
"""

import pandas as pd
from libsbml import Model as libModel
from refinegems.modelseed import get_modelseed_compounds

__author__ = "Famke Baeuerle"


def correct_charges_from_db(model: libModel, compounds: pd.DataFrame) -> tuple[libModel, dict]:
    """Adds charges taken from given database to metabolites which have no defined charge

    Args:
        - model (libModel): Model loaded with libsbml
        - compounds (pd.DataFrame): Containing database data with 'BiGG' (BiGG-Ids) and 'charge' (float or int) as columns

    Returns:
        tuple: libSBML model (1) & dictionary 'metabolite_id': list(charges) (2)
            (1) libModel: Model with added charges
            (2) dict: Metabolites with respective multiple charges
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


def correct_charges_modelseed(model: libModel) -> tuple[libModel, dict]:
    """Wrapper function which completes the steps to charge correction with the ModelSEED database

    Args:
        - model (libModel): Model loaded with libsbml

    Returns:
        tuple: libSBML model (1) & dictionary 'metabolite_id': list(charges) (2)
            (1) libModel: Model with added charges
            (2) dict: Metabolites with respective multiple charges
    """
    modelseed_compounds = get_modelseed_compounds()
    model_corr, multiple_charges = correct_charges_from_db(
        model, modelseed_compounds)
    
    return model_corr, multiple_charges
