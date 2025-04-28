"""Collection of utility functions."""

__author__ = "Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import bioregistry
import cobra
import logging

import memote.support.helpers as helpers
from memote.utils import truncate

from libsbml import Reaction
from six import iteritems
from typing import Union

################################################################################
# variables
################################################################################

# SBO terms
# ---------
SBO_BIOCHEM_TERMS = [
    "SBO:0000377",
    "SBO:0000399",
    "SBO:0000402",
    "SBO:0000403",
    "SBO:0000660",
    "SBO:0000178",
    "SBO:0000200",
    "SBO:0000214",
    "SBO:0000215",
    "SBO:0000217",
    "SBO:0000218",
    "SBO:0000219",
    "SBO:0000220",
    "SBO:0000222",
    "SBO:0000223",
    "SBO:0000233",
    "SBO:0000376",
    "SBO:0000401",
]  #: :meta:

SBO_TRANSPORT_TERMS = [
    "SBO:0000658",
    "SBO:0000657",
    "SBO:0000654",
    "SBO:0000659",
    "SBO:0000660",
]  #: :meta:

# Database local identifier regex patterns
# ----------------------------------------
DB2REGEX = bioregistry.get_pattern_map()  #: :meta hide-value:

# compartments
# ------------
VALID_COMPARTMENTS = {
    "c": "cytosol",
    "e": "extracellular space",
    "p": "periplasm",
    "uc": "unknown compartment",
}  #: :meta:
COMP_MAPPING = {
    # standard names
    "c": "c",
    "e": "e",
    "p": "p",
    # C_ naming
    "C_c": "c",
    "C_e": "e",
    "C_p": "p",
    # unknown compartments
    "": "uc", 
    "w": "uc",
}  #: :meta:

# useful defaults
# ---------------
MIN_GROWTH_THRESHOLD = 1.0e-5  #: :meta:

################################################################################
# functions
################################################################################

# SBO
# ---


def reannotate_sbo_memote(model: cobra.Model) -> cobra.Model:
    """Reannotate the SBO annotations (e.g. from SBOannotator) of a model
    into the SBO scheme accessible by memote.

    Args:
       - model (cobra.Model):
          The cobra Model to be reannotated.

    Returns:
        cobra.Model:
          The reannotated model
    """

    # reactions
    for r in model.reactions:
        if "sbo" in r.annotation:
            # biochemical
            if r.annotation["sbo"] in SBO_BIOCHEM_TERMS:
                r.annotation["sbo"] = "SBO:0000176"
            # transport
            elif r.annotation["sbo"] in SBO_TRANSPORT_TERMS:
                r.annotation["sbo"] = "SBO:0000185"
    # no change needed, as as of now, only one term exists
    # exchange
    # memote: SBO:0000627 / same for SBOannotator

    # demand
    # memote: SBO:0000628 / same for SBOannotator

    # sink
    # memote: SBO:0000632 / same for SBOannotator


# handling stoichiometric factors
# -------------------------------
def is_stoichiometric_factor(s: str) -> bool:
    """ "Check if a string could be used as a stoichiometric factor.

    Args:
        - s (str):
            The string to check.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


# handling biomass reaction
# -------------------------
def test_biomass_presence(model: cobra.Model) -> Union[list[str], None]:
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
    if outcome:
        return biomass_rxn
    else:
        return None


def sum_biomass_weight(reaction: Reaction) -> float:
    """
    From MEMOTE: https://github.com/opencobra/memote/blob/81a55a163262a0e06bfcb036d98e8e551edc3873/src/memote/support/biomass.py#L95

    Compute the sum of all reaction compounds.

    This function expects all metabolites of the biomass reaction to have
    formula information assigned.

    Args:
        - reaction (Reaction):
            The biomass reaction of the model under investigation.

    Returns:
        float:
            The molecular weight of the biomass reaction in units of g/mmol.
    """
    return (
        sum(
            -coef * met.formula_weight
            for (met, coef) in iteritems(reaction.metabolites)
        )
        / 1000.0
    )


def test_biomass_consistency(model: cobra.Model, reaction_id: str) -> Union[float, str]:
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

    Args:
        - model(cobraModel):
            The model loaded with COBRApy.
        - reaction_id(str):
            Reaction ID of a BOF.

    Returns:
        (1) Case: problematic input
                str:
                    an error message.

        (2) Case: successful testing
                float:
                    biomass weight
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
        if (1 - 1e-03) < biomass_weight < (1 + 1e-06):
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
    # outcome = (1 - 1e-03) < biomass_weight < (1 + 1e-06) -> Need to implement that for check
    # To account for numerical inaccuracies, a range from 1-1e0-3 to 1+1e-06
    # is implemented in the assertion check
    return biomass_weight
