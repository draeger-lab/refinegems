#!/usr/bin/env python
"""Most functions within this module were copied from the MEMOTE GitHub page and modified by Gwendolyn O. Döbel.

This module provides functions to be used to assess the biomass weight as well as normalise it.
"""

import os
import cobra
import logging
from cobra import Reaction
from cobra import Model as cobraModel
from refinegems.io import load_model_libsbml
from six import iteritems
import memote.support.helpers as helpers
from memote.utils import truncate
from typing import Union

__author__ = "MEMOTE and Gwendolyn O. Döbel"


def test_biomass_presence(model: cobraModel) -> Union[list[str], None]:
    """
    Modified from MEMOTE: https://github.com/opencobra/memote/blob/81a55a163262a0e06bfcb036d98e8e551edc3873/src/memote/suite/tests/test_biomass.py#LL42C3-L42C3
    
    Expect the model to contain at least one biomass reaction.

    The biomass composition aka biomass formulation aka biomass reaction
    is a common pseudo-reaction accounting for biomass synthesis in
    constraints-based modelling. It describes the stoichiometry of
    intracellular compounds that are required for cell growth. While this
    reaction may not be relevant to modeling the metabolism of higher
    organisms, it is essential for single-cell modeling.

    Implementation:
    Identifies possible biomass reactions using two principal steps:
    
        1. Return reactions that include the SBO annotation "SBO:0000629" for
        biomass.
        
        2. If no reactions can be identified this way:
        
            1. Look for the ``buzzwords`` "biomass", "growth" and "bof" in reaction IDs.
            2. Look for metabolite IDs or names that contain the ``buzzword`` "biomass" and obtain the set of reactions they are involved in.
            3. Remove boundary reactions from this set.
            4. Return the union of reactions that match the buzzwords and of the reactions that metabolites are involved in that match the buzzword.
        
    This test checks if at least one biomass reaction is present.
    
    If no reaction can be identified return None.

    """
    biomass_rxn = [rxn.id for rxn in helpers.find_biomass_reaction(model)]
    outcome = len(biomass_rxn) > 0
    logging.info(
        """In this model the following {} biomass reaction(s) were
        identified: {}""".format(
            len(biomass_rxn), truncate(biomass_rxn)
        )
    )
    if outcome: return biomass_rxn
    else: return None


def sum_biomass_weight(reaction: Reaction) -> float:
    """
    From MEMOTE: https://github.com/opencobra/memote/blob/81a55a163262a0e06bfcb036d98e8e551edc3873/src/memote/support/biomass.py#L95
    
    Compute the sum of all reaction compounds.

    This function expects all metabolites of the biomass reaction to have
    formula information assigned.

    Args:
        - reaction (Reaction): The biomass reaction of the model under investigation.

    Returns:
        float: The molecular weight of the biomass reaction in units of g/mmol.

    """
    return (
        sum(
            -coef * met.formula_weight
            for (met, coef) in iteritems(reaction.metabolites)
        )
        / 1000.0
    )
    
    
def test_biomass_consistency(model: cobraModel, reaction_id: str) -> Union[float, str]:
    """
    Modified from MEMOTE: https://github.com/opencobra/memote/blob/81a55a163262a0e06bfcb036d98e8e551edc3873/src/memote/suite/tests/test_biomass.py#L89
    
    Expect biomass components to sum up to 1 g[CDW].

    This test only yields sensible results if all biomass precursor
    metabolites have chemical formulas assigned to them.
    The molecular weight of the biomass reaction in metabolic models is
    defined to be equal to 1 g/mmol. Conforming to this is essential in order
    to be able to reliably calculate growth yields, to cross-compare models,
    and to obtain valid predictions when simulating microbial consortia. A
    deviation from 1 - 1E-03 to 1 + 1E-06 is accepted.

    Implementation:
    Multiplies the coefficient of each metabolite of the biomass reaction with
    its molecular weight calculated from the formula, then divides the overall
    sum of all the products by 1000.

    """
    reaction = model.reactions.get_by_id(reaction_id)
    try:
        biomass_weight = sum_biomass_weight(reaction)
    except TypeError:
        message = """
        One or more of the biomass components do not have a defined formula or contain unspecified chemical groups.
        The biomass overall weight could thus not be calculated.
        """
        return message
    else:
        if ((1 - 1e-03) < biomass_weight < (1 + 1e-06)):
            logging.info(            
                """The component molar mass of the biomass reaction {} sums up to {}
                which is inside the 1e-03 margin from 1 mmol / g[CDW] / h.
                """.format(
                    reaction_id, biomass_weight
                    )
                )
        else:
            logging.warning(
                """The component molar mass of the biomass reaction {} sums up to {}
                which is outside of the 1e-03 margin from 1 mmol / g[CDW] / h.
                """.format(
                    reaction_id, biomass_weight
                )
            )
    #outcome = (1 - 1e-03) < biomass_weight < (1 + 1e-06) -> Need to implement that for check
    # To account for numerical inaccuracies, a range from 1-1e0-3 to 1+1e-06
    # is implemented in the assertion check
    return biomass_weight


def normalise_biomass(biomass: Reaction, current_sum: float) -> Reaction:
    """Normalises the coefficients according to current biomass weight to one g[CDW]

    Args:
        - biomass (Reaction): Biomass function/reaction
        - current_sum (float): Biomass weight calculated with sum_biomass_weight in g/mmol

    Returns:
        Reaction: Biomass function/reaction with updated coefficients
    """
    metabs = biomass.metabolites # Get all metabolites
    
    # Normalise & update coefficients
    for (met, coef) in iteritems(metabs):
        metabs[met] = 1/current_sum * coef
    biomass._metabolites = metabs
    
    return biomass

def check_normalise_biomass(model: cobraModel) -> Union[cobraModel, None]:
    """
       1. Checks if at least one biomass reaction is present
       2. For each found biomass reaction checks if it sums up to 1g[CDW]
       3. Normalises the coefficients of each biomass reaction where the sum is not 1g[CDW] until the sum is 1g[CDW]
       4. Returns model with adjusted biomass function(s)
       

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        cobraModel: COBRApy model with adjusted biomass functions
    """
    
    biomass_rxn = test_biomass_presence(model)
    
    if biomass_rxn:
        for bm_rxn in biomass_rxn:
            bm_weight = test_biomass_consistency(model, bm_rxn)
            if type(bm_weight) == str: logging.error(f'Reaction {bm_rxn}: {bm_weight}')
            else:
                while not ((1 - 1e-03) < bm_weight < (1 + 1e-06)):
                    #normalise_biomass(model.reactions.get_by_id(bm_rxn), bm_weight)
                    model.reactions.get_by_id(bm_rxn).__imul__(1/bm_weight) # -> Maybe add like this?
                    bm_weight = test_biomass_consistency(model, bm_rxn)
                    logging.info(f'For reaction \'{bm_rxn}\' the coefficients changed.')
                
        cobra.io.write_sbml_model(model, f'../{model.id}_tmp.xml')
        model = load_model_libsbml(f'../{model.id}_tmp.xml')
        os.remove(f'../{model.id}_tmp.xml')
            
        return model
    
    else:
        logging.error(f'No biomass objective function was found in the provided model {model.id}.')
        return biomass_rxn
