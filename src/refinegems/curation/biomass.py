#!/usr/bin/env python
"""Most functions within this module were copied from the MEMOTE GitHub page and modified by Gwendolyn O. Döbel.

This module provides functions to normalise the biomass objective function(s).
"""

__author__ = "MEMOTE and Gwendolyn O. Döbel"

############################################################################
# requirements
############################################################################

import logging

from cobra import Model as cobraModel
from cobra import Reaction
from six import iteritems
from typing import Union

from ..utility.util import test_biomass_consistency, test_biomass_presence

############################################################################
# variables
############################################################################

############################################################################
# functions
############################################################################


# @DELETE - one the function below works?
def normalise_biomass(biomass: Reaction, current_sum: float) -> Reaction:
    """Normalises the coefficients according to current biomass weight to one g[CDW]

    Args:
        - biomass (Reaction): 
            Biomass function/reaction
        - current_sum (float): 
            Biomass weight calculated with sum_biomass_weight in g/mmol

    Returns:
        Reaction: 
            Biomass function/reaction with updated coefficients
    """
    metabs = biomass.metabolites # Get all metabolites
    
    # Normalise & update coefficients
    for (met, coef) in iteritems(metabs):
        metabs[met] = 1/current_sum * coef
    biomass._metabolites = metabs
    
    return biomass

def check_normalise_biomass(model: cobraModel, cycles:int=10) -> Union[cobraModel, None]:
    """
       1. Checks if at least one biomass reaction is present
       2. For each found biomass reaction checks if it sums up to 1g[CDW]
       3. Normalises the coefficients of each biomass reaction where the sum is not 1g[CDW] until the sum is 1g[CDW]
       4. Returns model with adjusted biomass function(s)
       

    Args:
        - model (cobraModel): 
            Model loaded with COBRApy
        - cycles (int, optional): 
            Maximal number of optiomisation cycles that will be run.
            Used to avoid endless optiomisation cycles.

    Returns:
        cobraModel: 
            COBRApy model with adjusted biomass functions
    """
    
    biomass_rxn = test_biomass_presence(model)
    
    if biomass_rxn:
        for bm_rxn in biomass_rxn:
            bm_weight = test_biomass_consistency(model, bm_rxn)
            if type(bm_weight) == str: logging.error(f'Reaction {bm_rxn}: {bm_weight}')
            else:
                c = 0 # counter to ensure it does not run endlessly
                while not ((1 - 1e-03) < bm_weight < (1 + 1e-06)) and c <= cycles:
                    #normalise_biomass(model.reactions.get_by_id(bm_rxn), bm_weight)
                    model.reactions.get_by_id(bm_rxn).__imul__(1/bm_weight) # -> Maybe add like this?
                    bm_weight = test_biomass_consistency(model, bm_rxn)
                    logging.info(f'For reaction \'{bm_rxn}\' the coefficients changed.')
                    c += 1
            
        return model
    
    else:
        logging.error(f'No biomass objective function was found in the provided model {model.id}.')
        return biomass_rxn
