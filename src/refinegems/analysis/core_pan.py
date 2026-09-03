"""Handling, creating and working with pan-core models."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra

from typing import Literal

from ..utility.io import load_model
from ..classes.reports import CorePanAnalysisReport
from ..utility.entities import resolve_compartment_names

################################################################################
# functions
################################################################################

# core-pan modelling
# ------------------


def extract_reactions_ids(
    model: cobra.Model, based_on: Literal["id"] = "id"
) -> list[str]:
    """Extract reactions identifiers from a model.

    Based on:

    - id: extracts the actual IDs as set in the model.

    Args:
        - model (cobra.Model):
            The model to extract the IDs from.
            Loaded with COBRApy.
        - based_on (Literal['id'], optional):
            How and which IDs to extract. Defaults to 'id'.

    Raises:
        - ValueError: Unknown input for parameter based_on if not in given options.

    Returns:
        list[str]:
            List of extracted IDs in the given format.
    """

    match based_on:
        case "id":
            return [_.id for _ in model.reactions]
        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")


def find_core_reaction_ids(all_reactions: dict[str : list[str]]) -> list[str]:
    """Helper function for :py:func:`~refinegems.analysis.core_pan.generate_core_pan_model`.
    Identify the core reactions from a set of reactions from different models.
    Core reactions are reactions that occur in ALL the models.

    Args:
        - all_reactions (dict[str:list[str]]):
            List of reactions IDs for all model to be part of the core-pan model.

    Returns:
        list[str]:
            List of the IDs of reactions that are defined as core.
    """

    core = []
    first = True
    for reacs in all_reactions.values():
        if first:
            core = set(reacs)
            first = False
        else:
            core = core.intersection(set(reacs))

    return core


def find_pan_reactions(
    all_reactions: dict[str : list[str]], core: list[str]
) -> list[str]:
    """Helper function for :py:func:`~refinegems.analysis.core_pan.generate_core_pan_model`. Identify the pan reactions
    for a set of reactions of different model. Pan reactions are reactions, that are found
    in AT LEAST one model but NOT in all.

    Args:
        - all_reactions (dict[str:list[str]):
            List of reactions IDs for all model to be part of the core-pan model.
        - core (list[str]):
            List of core reaction IDs, output of :py:func:`~refinegems.analysis.core_pan.find_core_reaction_ids`.

    Returns:
        list[str]:
            List of pan reaction IDs.
    """

    pan = {}
    for model, reacs in all_reactions.items():
        pan[model] = set([_ for _ in reacs if _ not in core])

    return pan


def collect_reacs_from_model(
    model: cobra.Model,
    reac_id_list: list[str],
    based_on: Literal["id"] = "id",
    notes: tuple[str] = ("core-pan", "core"),
) -> list[cobra.Reaction]:
    """Based on a model and a list of reactions IDs, collects the corresponding reactions.

    Args:
        - model (cobra.Model):
            The model.
        - reac_id_list (list[str]):
            List of reactions IDs. are treated as actual cobra ID or not depending on 'based_on'.
        - based_on (Literal['id'], optional):
            Defines, if the IDs are to be treated literal ('id') or not.
            Defaults to 'id'.
        - notes (tuple, optional):
            What kind of reactions have been collected. Expects a tuple of two strings.
            Uses the tuple to create a notes entry in the reaction object.
            Defaults to ('core-pan','core').

    Raises:
        - ValueError: Unknown input for parameter based_on.

    Returns:
        list[cobra.Reaction]:
            List of the extracted reactions.
    """

    match based_on:
        case "id":

            reac_list = []
            for id in reac_id_list:
                new_reac = model.reactions.get_by_id(id)
                new_reac.notes[notes[0]] = notes[1]
                reac_list.append(new_reac)

            return reac_list

        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")


def generate_core_pan_model(
    model_list: list[str],
    based_on: Literal["id"] = "id",
    name: str = "core_pan_model",
    remove_genes: bool = True,
) -> cobra.Model:
    """Generate a core-pan model from a set of models.

    Generation id based on:

        - id: uses the IDs to compare reactions

    Args:
        - model_list (list[str]):
            List of paths to models.
        - based_on (Literal['id'], optional):
            How to decide which reactions are considered the same.
            Defaults to 'id'.
        - name (str, optional):
            Name of the new model.
            Defaults to 'core_pan_model'.
        - remove_genes (bool, optional):
            Flag to remove all genes from the model.
            Defaults to True.

    Returns:
        cobra.Model:
            The generated core-pan model.
    """

    # load all models
    all_models = load_model(model_list, "cobra")

    # resolve compartment issue
    for model in all_models:
        resolve_compartment_names(model)

    # extract reactions
    all_reactions = {
        model.id: extract_reactions_ids(model, based_on) for model in all_models
    }

    # define core-pan
    core = find_core_reaction_ids(all_reactions)
    pan = find_pan_reactions(all_reactions, core)

    # extract corresponding reactions from input models
    core_reacs = collect_reacs_from_model(
        all_models[0], core, based_on, notes=("core-pan", "core")
    )
    pan_reacs = []
    collected = []
    for model, reacs in pan.items():
        to_add = [_ for _ in reacs if _ not in collected]
        current_model = [_ for _ in all_models if _.id == model][0]
        new_reacs = collect_reacs_from_model(
            current_model, to_add, based_on, notes=("core-pan", "pan")
        )
        pan_reacs.extend(new_reacs)
        collected.extend(to_add)

    # construct model
    cp_model = cobra.Model(name)
    cp_model.add_reactions(core_reacs)
    cp_model.add_reactions(pan_reacs)

    # step 4: remove genes (optional)
    if remove_genes:
        cobra.manipulation.delete.remove_genes(
            cp_model, cp_model.genes, remove_reactions=False
        )

    return cp_model


# core-pan comparison
# -------------------

def compare_to_core_pan(
    model: cobra.Model, cp_model: cobra.Model, based_on: Literal["id"] = "id"
) -> CorePanAnalysisReport:
    """Compare a model to a pan-core model.

    Comparison can be done based on:

        - id: uses the reaction IDs for a simple and direct comparison.
        
        .. note:: 
            Currently, this requires the model reactions to be annotated with 'core-pan' notes.
            This function however, is object to change and will be extended in the future.

    Args:
        - model (cobra.Model):
            The input model.
        - cp_model (cobra.Model):
            The core-pan model
        - based_on (Literal['id'], optional):
            How to perform the comparison.
            Defaults to 'id'.

    Raises:
        - ValueError: Unknown input for parameter based_on.

    Returns:
        CorePanAnalysisReport:
            The analysis results in form of a report object.
    """

    results = CorePanAnalysisReport(model)

    match based_on:
        # compare models solely based on their reactions IDs
        case "id":

            # separate cp_model reactions into core and pan list
            core_reac_list = [
                _.id for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "core"
            ]
            pan_reac_list = [
                _.id for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "pan"
            ]

            # compare model to the core and pan reaction list
            results.core_reac = [
                _.id for _ in model.reactions if _.id in core_reac_list
            ]
            results.pan_reac = [_.id for _ in model.reactions if _.id in pan_reac_list]
            results.novel_reac = [
                _.id
                for _ in model.reactions
                if _.id not in results.pan_reac and _.id not in results.core_reac
            ]

        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")

    return results
