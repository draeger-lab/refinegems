"""Handling, creating and working with pan-core models."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra

from typing import Literal, Optional

from ..utility.io import load_model
from ..classes.reports import CorePanAnalysisReport
from ..utility.entities import resolve_compartment_names

################################################################################
# functions
################################################################################

# core-pan modelling
# ------------------


def extract_reactions_ids(
    model: cobra.Model, based_on: Literal["id", "annotation"] = "id", db_prefix: Optional[str] = None
) -> list[str]:
    """Extract reactions identifiers from a model.

    Based on:

    - id: extracts the actual IDs as set in the model.
    - annotation: extracts the database ID from the reaction annotations based on db_prefix.

    Args:
        - model (cobra.Model):
            The model to extract the IDs from.
            Loaded with COBRApy.
        - based_on (Literal['id', 'annotation'], optional):
            How and which IDs to extract. Defaults to 'id'.
        - db_prefix (str, optional):
            The annotation database prefix to extract (e.g., 'kegg.reaction'). 
            Required if based_on is 'annotation'.

    Raises:
        - ValueError: Unknown input for parameter based_on if not in given options.
        - ValueError: Missing db_prefix when based_on is 'annotation'.

    Returns:
        list[str]:
            List of extracted IDs in the given format.
    """

    match based_on:
        case "id":
            return [_.id for _ in model.reactions]
        case "annotation":
            if not db_prefix:
                raise ValueError("db_prefix must be provided when based_on is 'annotation'")
            extracted = []
            for _ in model.reactions:
                ann = _.annotation.get(db_prefix)
                if isinstance(ann, list) and len(ann) > 0:
                    extracted.append(ann[0])
                elif isinstance(ann, str):
                    extracted.append(ann)
                else:
                    # Fallback to model ID if annotation is missing to prevent data loss
                    extracted.append(_.id) 
            return extracted
        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")


def find_core_reaction_ids(all_reactions: dict[str, list[str]]) -> list[str]:
    """Helper function for :py:func:`~refinegems.analysis.core_pan.generate_core_pan_model`.
    Identify the core reactions from a set of reactions from different models.
    Core reactions are reactions that occur in ALL the models.

    Args:
        - all_reactions (dict[str, list[str]]):
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

    return list(core)


def find_pan_reactions(
    all_reactions: dict[str, list[str]], core: list[str]
) -> dict[str, list[str]]:
    """Helper function for :py:func:`~refinegems.analysis.core_pan.generate_core_pan_model`. Identify the pan reactions
    for a set of reactions of different model. Pan reactions are reactions, that are found
    in AT LEAST one model but NOT in all.

    Args:
        - all_reactions (dict[str, list[str]]):
            List of reactions IDs for all model to be part of the core-pan model.
        - core (list[str]):
            List of core reaction IDs, output of :py:func:`~refinegems.analysis.core_pan.find_core_reaction_ids`.

    Returns:
        dict[str, list[str]]:
            Dictionary of pan reaction IDs per model.
    """

    pan = {}
    for model, reacs in all_reactions.items():
        pan[model] = list(set([_ for _ in reacs if _ not in core]))

    return pan


def collect_reacs_from_model(
    model: cobra.Model,
    reac_id_list: list[str],
    based_on: Literal["id", "annotation"] = "id",
    notes: tuple[str, str] = ("core-pan", "core"),
    db_prefix: Optional[str] = None
) -> list[cobra.Reaction]:
    """Based on a model and a list of reactions IDs, collects the corresponding reactions.

    Args:
        - model (cobra.Model):
            The model.
        - reac_id_list (list[str]):
            List of reactions IDs. are treated as actual cobra ID or not depending on 'based_on'.
        - based_on (Literal['id', 'annotation'], optional):
            Defines, if the IDs are to be treated literal ('id') or not.
            Defaults to 'id'.
        - notes (tuple, optional):
            What kind of reactions have been collected. Expects a tuple of two strings.
            Uses the tuple to create a notes entry in the reaction object.
            Defaults to ('core-pan','core').
        - db_prefix (str, optional):
            The annotation database prefix. Required if based_on is 'annotation'.

    Raises:
        - ValueError: Unknown input for parameter based_on.
        - ValueError: Missing db_prefix when based_on is 'annotation'.

    Returns:
        list[cobra.Reaction]:
            List of the extracted reactions.
    """

    reac_list = []
    match based_on:
        case "id":
            for id in reac_id_list:
                try:
                    new_reac = model.reactions.get_by_id(id)
                    new_reac.notes[notes[0]] = notes[1]
                    reac_list.append(new_reac)
                except KeyError:
                    pass

        case "annotation":
            if not db_prefix:
                raise ValueError("db_prefix must be provided when based_on is 'annotation'")
            
            # Create a lookup dictionary mapping the extracted identity to the reaction object
            ann_to_reac = {}
            for _ in model.reactions:
                ann = _.annotation.get(db_prefix)
                key = ann[0] if isinstance(ann, list) and len(ann) > 0 else (ann if isinstance(ann, str) else _.id)
                ann_to_reac[key] = _
            
            for id_val in reac_id_list:
                if id_val in ann_to_reac:
                    new_reac = ann_to_reac[id_val]
                    new_reac.notes[notes[0]] = notes[1]
                    reac_list.append(new_reac)

        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")
            
    return reac_list


def generate_core_pan_model(
    model_list: list[str],
    based_on: Literal["id", "annotation"] = "id",
    name: str = "core_pan_model",
    remove_genes: bool = True,
    db_prefix: Optional[str] = None
) -> cobra.Model:
    """Generate a core-pan model from a set of models.

    Generation id based on:

        - id: uses the IDs to compare reactions
        - annotation: uses database annotations to compare reactions

    Args:
        - model_list (list[str]):
            List of paths to models.
        - based_on (Literal['id', 'annotation'], optional):
            How to decide which reactions are considered the same.
            Defaults to 'id'.
        - name (str, optional):
            Name of the new model.
            Defaults to 'core_pan_model'.
        - remove_genes (bool, optional):
            Flag to remove all genes from the model.
            Defaults to True.
        - db_prefix (str, optional):
            The annotation database prefix. Required if based_on is 'annotation'.

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
        model.id: extract_reactions_ids(model, based_on, db_prefix) for model in all_models
    }

    # define core-pan
    core = find_core_reaction_ids(all_reactions)
    pan = find_pan_reactions(all_reactions, core)

    # extract corresponding reactions from input models
    core_reacs = collect_reacs_from_model(
        all_models[0], core, based_on, notes=("core-pan", "core"), db_prefix=db_prefix
    )
    
    pan_reacs = []
    collected = []
    for model, reacs in pan.items():
        to_add = [_ for _ in reacs if _ not in collected]
        current_model = [_ for _ in all_models if _.id == model][0]
        new_reacs = collect_reacs_from_model(
            current_model, to_add, based_on, notes=("core-pan", "pan"), db_prefix=db_prefix
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
    model: cobra.Model, 
    cp_model: cobra.Model, 
    based_on: Literal["id", "annotation"] = "id",
    db_prefix: Optional[str] = None
) -> CorePanAnalysisReport:
    """Compare a model to a pan-core model.

    Comparison can be done based on:

        - id: uses the reaction IDs for a simple and direct comparison.
        - annotation: uses database annotations for comparison.
        
        .. note:: 
            Currently, this requires the model reactions to be annotated with 'core-pan' notes.
            This function however, is object to change and will be extended in the future.

    Args:
        - model (cobra.Model):
            The input model.
        - cp_model (cobra.Model):
            The core-pan model
        - based_on (Literal['id', 'annotation'], optional):
            How to perform the comparison.
            Defaults to 'id'.
        - db_prefix (str, optional):
            The annotation database prefix. Required if based_on is 'annotation'.

    Raises:
        - ValueError: Unknown input for parameter based_on.
        - ValueError: Missing db_prefix when based_on is 'annotation'.

    Returns:
        CorePanAnalysisReport:
            The analysis results in form of a report object.
    """

    results = CorePanAnalysisReport(model)

    match based_on:
        case "id":
            core_reac_list = [
                _.id for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "core"
            ]
            pan_reac_list = [
                _.id for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "pan"
            ]

            results.core_reac = [_.id for _ in model.reactions if _.id in core_reac_list]
            results.pan_reac = [_.id for _ in model.reactions if _.id in pan_reac_list]
            results.novel_reac = [
                _.id for _ in model.reactions if _.id not in results.pan_reac and _.id not in results.core_reac
            ]

        case "annotation":
            if not db_prefix:
                raise ValueError("db_prefix must be provided when based_on is 'annotation'")
            
            def get_identity(_reac):
                ann = _reac.annotation.get(db_prefix)
                return ann[0] if isinstance(ann, list) and len(ann) > 0 else (ann if isinstance(ann, str) else _reac.id)

            core_reac_list = [
                get_identity(_) for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "core"
            ]
            pan_reac_list = [
                get_identity(_) for _ in cp_model.reactions if "core-pan" in _.notes and _.notes["core-pan"] == "pan"
            ]

            results.core_reac = []
            results.pan_reac = []
            results.novel_reac = []

            for _ in model.reactions:
                identity = get_identity(_)
                if identity in core_reac_list:
                    results.core_reac.append(_.id)
                elif identity in pan_reac_list:
                    results.pan_reac.append(_.id)
                else:
                    results.novel_reac.append(_.id)

        case _:
            raise ValueError(f"Unknown input for parameter based_on: {based_on}")

    return results