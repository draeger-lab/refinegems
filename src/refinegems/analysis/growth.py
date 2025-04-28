#!/usr/bin/env python
"""Provides functions to simulate growth on any medium and other functionalities replated to growth."""

__author__ = "Famke Baeuerle and Carolin Brune"

############################################################################
# requirements
############################################################################

import cobra
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

from cobra import Model as cobraModel
from ..utility.util import test_biomass_presence
from ..utility.io import load_model, load_a_table_from_database
from ..classes.reports import (
    SingleGrowthSimulationReport,
    GrowthSimulationReport,
    AuxotrophySimulationReport,
    SourceTestReport,
)
from ..classes.medium import (
    Medium,
    medium_to_model,
    read_from_cobra_model,
    load_medium_from_db,
    load_media,
)
from typing import Literal, Union

############################################################################
# functions
############################################################################

# growth simulation and more
# --------------------------


def set_bounds_to_default(
    model: cobraModel, reac_bounds: Union[None, str, tuple[float]] = None
):
    """Set the reactions bounds of a model to given default values.
    (Ir)reversibility is retained.

    Args:
        - model (cobraModel):
            The model loaded with COBRApy.
        - reac_bounds (None|str|tuple[float], optional):
            The setting for the new reaction bounds.
            Defaults to None. If None or "cobra", uses the COBRApy in-built default values (-1000.0, 1000.0).
            The user can set personal values by entering a tuple of two floats.

    Raises:
        - ValueError: Problematic input for bounds, if neither None, "cobra" or a tuple of floats is entered for reac_bounds.
    """

    # user-specific default bounds (tuple of two floats)
    if (
        type(reac_bounds) is tuple
        and len(reac_bounds) == 2
        and type(reac_bounds[0]) is float
        and type(reac_bounds[1]) is float
    ):
        pass
    # use COBRApy-internal default bounds (None or string 'cobra')
    elif reac_bounds is None or (type(reac_bounds) is str and reac_bounds == "cobra"):
        reac_bounds = cobra.Configuration().bounds
    # unknown input
    else:
        raise ValueError(f"Problematic input for bounds: {reac_bounds}")

    # apply the bounds to the model
    for reaction in model.reactions:

        # reactions is originally disabled
        if reaction.upper_bound == 0.0 and reaction.lower_bound == 0.0:
            pass
        # forward only
        elif reaction.lower_bound == 0.0:
            reaction.upper_bound = reac_bounds[1]
        # backward only
        elif reaction.upper_bound == 0.0:
            reaction.lower_bound = reac_bounds[0]
        # reversible or broken
        else:
            reaction.bounds = reac_bounds


def get_uptake(model: cobraModel, type: str) -> list[str]:
    """Compute the list of exchange reactions that have fluxes > 0 under certain conditions.

    Args:
        - model (cobraModel):
            A cobra Model to be tested.
        - type (str):
            Type of uptake, can be 'minimal'/'min' or 'standard'/'std'.

    Raises:
        - ValueError: Unknown type for uptake, if type not in ['minimal','min','standard','std']

    Returns:
        list[str]:
            List of non-zero flux exchange reactions under the set type.
    """

    match type:
        # return minimal
        case "minimal" | "min":
            with model:
                minimal = cobra.medium.minimal_medium(model)
                # print(minimal)
                return list(minimal.index)
        # return standart, non-zero flux compounds
        case "standard" | "std":
            with model:
                sol = model.optimize()
                fluxes = sol.fluxes
                uptake = []
                for index, value in fluxes.items():
                    if index in model.exchanges:
                        if value < 0:
                            uptake.append(index)
            return uptake
        case _:
            raise ValueError(f"Unknown type for uptake: {type}")


def get_secretion(model: cobraModel) -> list[str]:
    """Returns the list of exchange reactions for compounds that are secreted in the current version of the model.

    Args:
        - model (cobraModel):
            The cobra model to be tested.

    Returns:
        list[str]:
            The list of IDs of secretion reactions
    """

    with model:

        sf = model.summary().secretion_flux
        s = sf[sf["flux"] < 0.0].index.tolist()

    return s


def get_production(model: cobraModel) -> list[str]:
    """Checks fluxes after FBA, if positive the metabolite is produced.

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        list[str]:
            Ids of produced metabolites
    """

    with model:
        fluxes = model.optimize().fluxes
        production = []
        for index, value in fluxes.items():
            if value > 0:
                production.append(index)
    return production


def find_growth_essential_exchanges(
    model: cobraModel, growth_medium: dict, standard_uptake: Union[list[str], None]
) -> list[str]:
    """Find exchanges in a medium (with or without supplements) essential for the growth.

    .. note::

        This function currently only tests single deletions.

    Args:
        - model (cobraModel):
            The model to be tested.
        - growth_medium (dict):
            The medium in dictionary form, ready to be added to the model.
        - standard_uptake (list[str]|None):
            Option to add a second medium list as supplements.

    Returns:
        list[str]:
            The list of exchanges essential for growth.
    """
    with model:
        if standard_uptake:
            # combine standard and second medium, set all fluxes to 10.0 for standard
            standard_medium = {i: 10.0 for i in standard_uptake}
            new_medium = {**growth_medium, **standard_medium}
        else:
            new_medium = growth_medium
        # add new medium to model
        try:
            model.medium = new_medium
        except ValueError:
            logging.info(
                "Change upper bounds to COBRApy defaults to make model simulatable."
            )
            set_bounds_to_default(model)
            model.medium = new_medium
        # find essentials
        essential = []
        for metab in new_medium.keys():
            with model:
                model.reactions.get_by_id(metab).lower_bound = 0
                sol = model.optimize()
                if (
                    sol.objective_value < 1e-5
                ):  # and sol.objective_value > -1e-9: # == 0 no negative growth!
                    essential.append(metab)

    return essential


def find_additives_to_enable_growth(
    model: cobraModel,
    growth_medium: dict,
    standard_uptake: list[str],
    combine: bool = False,
) -> Union[list[str], dict]:
    """Based on a new medium for growth and a standard one the model already growths on, find additives from the standard,
    which can be added to the new one to enable growths.

    .. warning::
        This function currently only checks via single deletions not via correlated ones.
        Might lead to problems in complecated cases, e.g. missing carbon source but
        multiple ones in the standard medium.

    Args:
        - model (cobraModel):
            The model to add the medium to.
        - growth_medium (dict):
            A medium definition, ready to be added to a COBRApy model.
            This is for the new medium for growth testing.
        - standard_uptake (list[str]):
            A list of exchange reactions, e.g. output of get_uptake.
            This is for the old medium where the model can grow on.
        - combine (bool, optional):
            Flag to directly combine the additives with the new medium or just return the additive's IDs.
            Defaults to False (returns reaction IDs).

    Returns:
        (1) Case ``combine = False``:
                list:
                    List of the exchange reaction IDs of the additives.

        (2) Case ``combine = True``:
                Medium:
                    Supplemented medium as a dictionary.
    """

    # find essential exchange reactions
    essential = find_growth_essential_exchanges(model, growth_medium, standard_uptake)

    # find the essential compounds not in the growth medium
    additives = [_ for _ in essential if _ not in growth_medium.keys()]

    # return ...
    if combine:
        # ... the supplemented medium
        supplements = {_: 10.0 for _ in additives}
        suppl_medium = {**growth_medium, **supplements}
        return suppl_medium
    else:
        # ... the list of suplements
        return additives


def get_metabs_essential_for_growth_wrapper(
    model: cobraModel, media: list[Medium], only_additives: bool = True
) -> dict:
    """
    Returns metabolites necessary for growth and not in media

    Args:
        - model (cobraModel):
            Model loaded with COBRApy
        - media (list[Medium]):
            Containing all media for which the growth essential metabolites not contained in the media should be returned
        - only_additives(bool, optional):
            Flag to only return the supplemented exchanges (True) or all essential ones (False).
            Defaults to True.

    Returns:
        dict:
            Information on different media which metabs are missing (key: name of medium, values: list of exchanges).
    """

    default_uptake = get_uptake(model, "std")
    ess = {}
    for medium in media:

        # convert to a cobrapy medium
        exported_m = medium_to_model(medium=medium, model=model, add=False)
        # get essentials
        essential = find_growth_essential_exchanges(model, exported_m, default_uptake)
        # check for only_additives flag
        if only_additives:
            essential = [_ for _ in essential if _ not in exported_m.keys()]
        # add to dict
        ess[medium.name] = essential

    return ess


def growth_sim_single(
    model: cobraModel,
    m: Medium,
    namespace: Literal["BiGG", "Name"] = "BiGG",
    supplement: Literal[None, "std", "min"] = None,
) -> SingleGrowthSimulationReport:
    """Simulate the growth of a model on a given medium.

    Args:
        - model (cobraModel):
            The model.
        - m (Medium):
            The medium.
        - supplement (Literal[None,'std','min'], optional):
            Flag to add additvites to the model to ensure growth. Defaults to None (no supplements).
            Further options include 'std' for standard uptake and 'min' for minimal uptake supplementation.

    Returns:
        SingleGrowthSimulationReport:
            Object with the simulation results
    """

    with model:

        # convert to a cobrapy medium
        exported_m = medium_to_model(
            medium=m,
            model=model,
            namespace=namespace,
            default_flux=10.0,
            replace=False,
            double_o2=False,
            add=False,
        )

        # supplement, if tag is set
        match supplement:
            case "std":
                uptake = get_uptake(model, "std")
                new_m = find_additives_to_enable_growth(
                    model, exported_m, uptake, combine=True
                )
            case "min":
                uptake = get_uptake(model, "min")
                new_m = find_additives_to_enable_growth(
                    model, exported_m, uptake, combine=True
                )
            case _:
                new_m = exported_m

        # add medium to model
        try:
            model.medium = new_m
        except ValueError:
            logging.info(
                "Change upper bounds to 1000.0 and lower bounds to -1000.0 to make model simulatable."
            )
            set_bounds_to_default(model)
            model.medium = new_m

        # simulate growth
        report = SingleGrowthSimulationReport(model_name=model.id, medium_name=m.name)
        report.growth_value = model.optimize().objective_value
        report.doubling_time = (
            (np.log(2) / report.growth_value) * 60 if report.growth_value != 0 else 0
        )
        report.additives = [_ for _ in new_m if _ not in exported_m]
        report.no_exchange = [
            _
            for _ in m.export_to_cobra(namespace=namespace).keys()
            if _ not in exported_m
        ]

    return report


def growth_sim_multi(
    models: Union[cobraModel, list[cobraModel]],
    media: Union[Medium, list[Medium]],
    namespace: Literal["BiGG", "Name"] = "BiGG",
    supplement_modes: Union[
        list[Literal["None", "min", "std"]], None, Literal["None", "min", "std"]
    ] = None,
) -> GrowthSimulationReport:
    """Simulate the growth of (at least one) models on (at least one) media.

    Args:
        - models (cobraModel | list[cobraModel]):
            A COBRApy model or a list of multiple.
        - media (Medium | list[Medium]):
            A refinegems Medium object or a list of multiple.
        - supplement_modes (list[Literal[None,'min','std']] | None | Literal[None, 'min', 'std'], optional):
            Option to supplement the media to enable growth.
            Default to None. Further options include a list with one entry for each medium or a string to set the same default for all.
            The string can be 'min', 'std' or None.

    Returns:
        GrowthSimulationReport:
            The compiled information of the simulation results.
    """

    # check input
    if type(models) != list:
        models = [models]
    if type(media) != list:
        media = [media]
    if type(supplement_modes) != list:
        supplement_modes = [supplement_modes] * len(media)

    # simulate the growth of the models on the different media
    report = GrowthSimulationReport()
    for mod in models:
        for med, supp in zip(media, supplement_modes):
            r = growth_sim_single(mod, med, namespace=namespace, supplement=supp)
            report.add_sim_results(r)

    return report


def growth_analysis(
    models: Union[cobra.Model, str, list[str], list[cobra.Model]],
    media: Union[Medium, list[Medium], str],
    namespace: Literal["BiGG"] = "BiGG",
    supplements: Union[
        None, list[Literal[None, "std", "min"]], Literal[None, "std", "min"]
    ] = None,
    retrieve: Literal["report", "plot", "both"] = "plot",
) -> Union[GrowthSimulationReport, plt.Figure, tuple]:
    """Perform a growth analysis

    Args:
        - models (cobra.Model | str | list[str] | list[cobra.Model]):
            Model(s) to be tested.
            Can be a COBRA model, a path to one, or a list of either type.
        - media (Medium | list[Medium] | str):
            Medium or media to be tested.
            Can be single or a list of medium objects or a path to a medium config file.
        - namespace (Literal['BiGG'], optional):
            Namespace of the model.
            Defaults to 'BiGG'.
        - supplements (None | list[Literal[None,'std','min']] | Literal[None,'std','min'], optional):
            Option to supplement media to enable growth. Can be None, std or min. If a single
            string is given, uses it for all media. Alternatively, a list with one entry per medium
            can be given to choose different options for different media.
            Defaults to None.
        - retrieve (Literal['report','plot','both'], optional):
            Determine, what to return.
            Options are the report object, the plotted graph or both (As a tuple).
            Defaults to 'plot'.

    Raises:
        - TypeError: Unknown or mixed types in model list.
        - KeyError: Empty list for models detected.
        - ValueError: Unknown input type for models.
        - TypeError: Unknown type found in media, should be list fo Medium.
        - ValueError: Unknown input for media.
        - ValueError: Unknown input for retrieve

    Returns:
        (1) Case: ``retrieve = report``:
                GrowthSimulationReport:
                    The generated report object.

        (2) Case: ``retrieve = plot``:
                plt.Figure:
                    The finished plot

        (3) Case: ``retrieve = both``:
                tuple:
                    The (0) report and the graphic (1).
    """

    # read-in all models into list
    # ----------------------------
    mod_list = []
    match models:
        case list():
            # list as input
            if len(models) > 0:
                # if list entries are paths
                if all(isinstance(_, str) for _ in models):
                    mod_list = load_model(models, package="cobra")
                # if list entries are already cobra.Models
                elif all(isinstance(_, cobra.Model) for _ in models):
                    mod_list = models
                else:
                    raise TypeError("Unknown or mixed types in model list.")
            else:
                raise KeyError("Empty list for models detected.")
        # single model as input
        case cobra.Model():
            mod_list = [models]
        # single string as input
        case str():
            mod = load_model(models, "cobra")
            mod_list = [mod]
        # unknown input
        case _:
            raise ValueError(f"Unknown input type for models: {type(models)}")

    # collect all media into list
    # ---------------------------
    match media:
        # single medium
        case Medium():
            media = [media]
        # list of media
        case list():
            if all(isinstance(_, Medium) for _ in media):
                pass
            else:
                raise TypeError(
                    "Unknown type found in media, should be list fo Medium."
                )
        # string - connection to YAML config file
        case str():
            media, supplements = load_media(media)
        # unknown input
        case _:
            raise ValueError(f"Unknown input for media: {media}")

    # run simulation
    # --------------
    report = growth_sim_multi(mod_list, media, namespace, supplements)

    # save / visualise report
    # -----------------------
    match retrieve:
        case "report":
            return report
        case "plot":
            return report.plot_growth()
        case "both":
            return (report, report.plot_growth())
        case _:
            raise ValueError(f"Unknown input for retrieve: {retrieve}")


def get_essential_reactions_via_single_knockout(model: cobraModel) -> list[str]:
    """Knocks out each reaction, if no growth is detected the reaction is seen as essential

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        list[str]:
            Ids of essential reactions
    """
    ess = []
    for reaction in model.reactions:
        with model as model:
            reaction.knock_out()
            sol = model.optimize().objective_value
            if sol <= 11:
                print(
                    "%s blocked (bounds: %s), new growth rate %f $"
                    % (reaction.id, str(reaction.bounds), sol)
                )
                ess.append(reaction.id)

    return ess


def get_essential_exchanges_via_bounds(model: cobraModel) -> list[str]:
    """Knocks out reactions by setting their bounds to 0, if no growth is detected the reaction is seen as essential


    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        list[str]:
            Ids of essential reactions
    """
    medium = model.medium
    ess = []
    for content in medium.keys():
        model.reactions.get_by_id(content).lower_bound = 0.0
        solution = model.optimize().objective_value
        if solution < 1e-9:
            ess.append(content)
        model.reactions.get_by_id(content).lower_bound = -10.0

    return ess


def find_growth_enhancing_exchanges(
    model: cobraModel, base_medium: dict
) -> pd.DataFrame:
    """Iterates through all exchanges to find metabolites that lead to a higher growth rate compared to the growth rate yielded on the base_medium

    Args:
        - model (cobraModel):
            Model loaded with COBRApy
        - base_medium (dict):
            Exchanges as keys and their flux bound as value (f.ex {'EX_glc__D_e' : 10.0})

    Returns:
        pd.DataFrame:
            Exchanges sorted from highest to lowest growth rate improvement
    """
    with model:
        model.medium = base_medium
        sol = model.optimize()
    base_growth = sol.objective_value
    print(base_growth)

    enhancement = {}
    for ex in list(model.exchanges):
        with model:
            base_medium[ex.id] = 10.0
            model.medium = base_medium
            sol = model.optimize()
            if sol.objective_value > base_growth:
                enhancement[ex.id] = sol.objective_value - base_growth
            base_medium.pop(ex.id)

    adds = pd.DataFrame(enhancement.items(), columns=["exchange", "diff"]).sort_values(
        by=["diff"], ascending=False
    )

    return adds


# auxotrophy simulation
# ---------------------


def test_auxotrophies(
    model: cobraModel,
    media_list: list[Medium],
    supplement_list: list[Literal[None, "min", "std"]],
    namespace: Literal["BiGG", "Name"] = "BiGG",
) -> AuxotrophySimulationReport:
    """Test for amino acid auxothrophies for a model and a list of media.

    Tests, if the model growths on the media and if and with what fluxes the
    20 proteinogenic amino acids are produced by temporarily adding a
    sink reaction for each of the amino acids to the model as the objective function.

    Args:
        - model (cobraModel):
            The model to be tested. Loaded with COBRApy.
        - media_list (list[Medium]):
            List of media to be tested.
        - supplement_list (list[Literal[None,'min','std']]):
            List of supplement modes for the media.
        - namespace (Literal['BiGG','Name'], optional):
            String for the namespace to be used for the model.
            Current options include 'BiGG', 'Name'.
            Defaults to 'BiGG'.

    Raises:
        - ValueError: Unknown input for namespace parameter.

    Returns:
        AuxotrophySimulationReport:
            The report for the test containing a table of the amino acids and the media names containing the simualted flux values.
    """

    results = {}

    # get amino acids from database
    amino_acids = Medium("aa").add_subset(subset_name="protAA")
    aa_list = set(amino_acids.substance_table["name"])

    # iterate over all media
    for med, supp in zip(media_list, supplement_list):
        auxotrophies = {}
        # then iterate over all amino acids
        for a in aa_list:
            entry = amino_acids.substance_table.loc[
                (amino_acids.substance_table["name"] == a)
                & (amino_acids.substance_table["db_type"].str.contains(namespace))
            ]

            with model as m:

                # export the medium
                exported_m = medium_to_model(
                    m,
                    med,
                    namespace=namespace,
                    default_flux=10.0,
                    replace=False,
                    double_o2=False,
                    add=False,
                )
                # supplement, if tag is set
                match supp:
                    case "std":
                        uptake = get_uptake(model, "std")
                        new_m = find_additives_to_enable_growth(
                            model, exported_m, uptake, combine=True
                        )
                    case "min":
                        uptake = get_uptake(model, "min")
                        new_m = find_additives_to_enable_growth(
                            model, exported_m, uptake, combine=True
                        )
                    case _:
                        new_m = exported_m
                # add medium to model
                m.medium = new_m

                # check namespace availability
                if len(entry) == 0:
                    warn_str = f"Amino acid {a} has no identifier for your chosen namespace {namespace}. Please contact support if you want to add one."
                    warnings.warn(warn_str)
                    growth_res = m.optimize()
                    auxotrophies[a] = (
                        growth_res.objective_value
                        if growth_res.status == "optimal"
                        else 0.0
                    )
                else:

                    # create and check IDs for the chosen namespace
                    internal_meta = ""
                    match namespace:
                        case "BiGG":
                            for np_id in entry["db_id"]:
                                if np_id + "_c" in [_.id for _ in m.metabolites]:
                                    internal_meta = np_id + "_c"
                                    break

                            exchange_reac = "EX_" + internal_meta
                            sink_reac = f"sink_{internal_meta}_tmp"
                        case _:
                            raise ValueError(
                                "Unknown namespace: {namespace}. Cannot create IDs."
                            )

                    # create a pseudo reaction -> a sink reaction for the amino acid
                    # to use as the new objective
                    if internal_meta == "":
                        warnings.warn(f"No identifier matched in cytosol for {a}.")
                    else:
                        m.add_boundary(
                            m.metabolites.get_by_id(internal_meta),
                            type="sink",
                            reaction_id=sink_reac,
                        )
                        m.objective = sink_reac
                        # if existent, close the exchange reaction
                        if exchange_reac in [_.id for _ in m.exchanges]:
                            m.reactions.get_by_id(exchange_reac).lower_bound = 0.0
                            m.reactions.get_by_id(exchange_reac).upper_bound = 0.0
                    # and calculate the new objective
                    growth_res = m.optimize()
                    auxotrophies[a] = (
                        growth_res.objective_value
                        if growth_res.status == "optimal"
                        else 0.0
                    )

            # add the current test results to the list of all results
            results[med.name] = auxotrophies

    report = AuxotrophySimulationReport(pd.DataFrame.from_dict(results))

    return report


# source test
# -----------


def test_growth_with_source(
    model: cobra.Model,
    element: str,
    substances: Union[None, str, list[str]] = None,
    medium: Union[None, str, Medium] = None,
    namespace: Literal["BiGG"] = "BiGG",
) -> SourceTestReport:
    """Test the growth of a model when switching out the source of a given chemical element for
    a set medium.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - element (str):
            The chemical symbol e.g., N for nitrogen, to change the sources for.
        - substances (None | str | list[str], optional):
            Substances to switch out in the medium.
            Can be a list of substance names present in the database, a subset name to be
            loaded from the database or None, which results in all substances in the database,
            that contain the element being tested as a source. Option None can potentially run
            a while.
            Defaults to None.
        - medium (None | str | Medium, optional):
            The medium to start with.
            The chosen medium ideally should have all other necessary elements needed for the model
            to grow.
            Defaults to None.
        - namespace (Literal['BiGG'], optional):
            The namespace to work on.
            Defaults to 'BiGG'.

    Raises:
        - KeyError: No growth function in model. Please add one beforehand.

    Returns:
        SourceTestReport:
            A report object with the results.
    """

    # validate input
    # model is required to have a growth function
    growth_funcs = test_biomass_presence(model)
    if growth_funcs:
        if any(_ in str(model.objective.expression) for _ in growth_funcs):
            pass
        else:
            warnings.warn(
                f"No growth functions set as objective, but growth function(s) detected. Setting objective to {growth_funcs[0]}"
            )
            model.objective = growth_funcs[0]
    else:
        raise KeyError("No growth function in model. Please add one beforehand.")

    # get the starting medium
    match medium:
        case str():
            current_medium = load_medium_from_db(medium)
        case Medium():
            current_medium = medium
        case _:
            current_medium = read_from_cobra_model(model)

    # save old medium settings
    origin_medium = model.medium

    # get the sources
    match substances:

        # case 1:
        # take a given subset from the database
        case str():
            temp_medium = Medium(
                "temp",
                pd.DataFrame(
                    columns=["name", "formula", "flux", "source", "db_id", "db_type"]
                ),
            )
            temp_medium.add_subset(substances)
            source_list = temp_medium.substance_table["name"].to_list()

        # case 2:
        # use the user given list - account for errors
        case list():
            source_list = substances

        # case 3:
        # download all possible options - may take some time
        case _:
            # get complete table
            substances = load_a_table_from_database("substance", query=False)
            # regex to find substances with element - element letter code NOT followed by ANY small letter
            element_regex = element + r"(?![a-z])"
            substances_mask = substances["formula"].str.contains(element_regex)
            substances_mask.fillna(value=False, inplace=True)
            substances = substances[substances_mask]
            source_list = substances.name.to_list()

    # perform the test
    results = []
    for s in source_list:
        # set new source
        current_medium.set_source(element, s)
        # add medium to model
        medium_to_model(model, current_medium, namespace=namespace, add=True)
        # simulate growth
        growth_value = model.optimize().objective_value
        results.append({"substance": s, "growth value": growth_value})

    results = SourceTestReport(pd.DataFrame(results), element, model.id)

    # set original medium for model
    model.medium = origin_medium

    return results


# minimal medium
# --------------


def model_minimal_medium(
    model: cobraModel,
    objective: Literal["flux", "medium", "exchanges"] = "flux",
    growth_rate: float = 0.5,
    open_exchanges: bool = False,
) -> Medium:
    """Get the minimal medium based on different objectives:

    - 'flux':      find the minimal fluxes based in current medium.
    - 'medium':    find the minimal number of compounds for the current medium.
    - 'exchanges': find the minimal number of compounds in a medium based on all avaiblae exchange reactions in the model.

    .. note::

        There may be multiple solutions for the minimisation, but only 1 will be returned.

    Args:
        - model (cobraModel):
            Model with a medium, that should be minimised.
        - objective (Literal[flux,medium,exchanges], optional):
            Objective for the minimisation task.
            Options listed above. Defaults to 'flux'.
        - growth_rate (float, optional):
            Minimum growth rate the model has to archieve.
            Defaults to 0.5. Only needed for objectives medium and exchanges.
        - open_exchanges (bool, optional):
            If set to True assigns large upper bound to all import reactions.
            Defaults to False.

            .. warning::

                Running `open_exchanges` on `True` can lead to infeasible runtimes.

    Raises:
        - ValueError: unknown objective.

    Returns:
        Medium:
            The medium that is a solution for the minimisation task.
    """

    # minimise the fluxes of the current medium
    if objective == "flux":
        max_growth = model.slim_optimize()
        min_medium = dict(cobra.medium.minimal_medium(model, max_growth))

    # minimise components of current medium
    elif objective == "medium":
        warnings.warn(
            "Warning: cobrapy.minimal_medium uses MIP formulation. This may take some time."
        )
        min_medium = dict(
            cobra.medium.minimal_medium(
                model,
                growth_rate,
                minimize_components=True,
                open_exchanges=open_exchanges,
            )
        )

    # get minimal number of medium components
    # based on exchange reaction possible in the model
    # note 1: can be time consuming
    # note 2: can lead to different results if run only once each time
    elif objective == "exchanges":
        # create cobra medium from all available exchange reactions
        ex_medium = {_.id: 1000.0 for _ in model.exchanges}
        # perform minimisation
        model.medium = ex_medium
        warnings.warn(
            "Warning: cobrapy.minimal_medium uses MIP formulation. This may take some time."
        )
        min_medium = dict(
            cobra.medium.minimal_medium(
                model,
                growth_rate,
                minimize_components=True,
                open_exchanges=open_exchanges,
            )
        )

    else:
        raise ValueError("Unknown objective for minimisation.")

    # create a Medium object from the minimal medium
    with model as tmp_model:
        tmp_model.medium = min_medium
        medium = read_from_cobra_model(tmp_model)

    return medium
