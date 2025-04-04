"""Identify, report and solve energy generating cycles (EGCs)."""

__author__ = "Carolin Brune and Tobias Fehrenbach"

################################################################################
# requirements
################################################################################

from ..analysis.growth import set_bounds_to_default
from .medium import Medium
from ..utility.util import test_biomass_presence, MIN_GROWTH_THRESHOLD

import cobra
import pandas as pd
import warnings

from functools import partial
from multiprocess import Pool
from tqdm import tqdm
from typing import Literal, Union

################################################################################
# variables
################################################################################

DISSIPATION_RXNS = {
    "ATP": {
        "ATP [Adenosine triphosphate]": -1,
        "Water [H2O]": -1,
        "ADP [Adenosine diphosphate]": 1,
        "Hydrogen [H(+)]": 1,
        "Phosphate [PO4(3-)]": 1,
    },
    "CTP": {
        "CTP [Cytidine triphosphate]": -1,
        "Water [H2O]": -1,
        "CDP [Cytidine diphosphate]": 1,
        "Hydrogen [H(+)]": 1,
        "Phosphate [PO4(3-)]": 1,
    },
    "GTP": {
        "GTP [Guanosine triphosphate]": -1,
        "Water [H2O]": -1,
        "GDP [Guanosine diphosphate]": 1,
        "Hydrogen [H(+)]": 1,
        "Phosphate [PO4(3-)]": 1,
    },
    "UTP": {
        "UTP [Uridine triphosphate]": -1,
        "Water [H2O]": -1,
        "UDP [Uridine diphosphate]": 1,
        "Hydrogen [H(+)]": 1,
        "Phosphate [PO4(3-)]": 1,
    },
    "ITP": {
        "ITP [Inosine triphosphate]": -1,
        "Water [H2O]": -1,
        "IDP [Inosine diphosphate]": 1,
        "Hydrogen [H(+)]": 1,
        "Phosphate [PO4(3-)]": 1,
    },
    "NADH": {
        "NADH [reduced Nicotinamide adenine dinucleotide]": -1,
        "Hydrogen [H(+)]": 1,
        "NAD [oxidized Nicotinamide adenine dinucleotide]": 1,
    },
    "NADPH": {
        "NADPH [reduced Nicotinamide adenine dinucleotide phosphate]": -1,
        "Hydrogen [H(+)]": 1,
        "NADP [oxidized Nicotinamide adenine dinucleotide phosphate]": 1,
    },
    "FADH2": {
        "FADH2 [reduced Flavin adenine dinucleotide]": -1,
        "Hydrogen [H(+)]": 2,
        "FAD [oxidized Flavin adenine dinucleotide]": 1,
    },
    "FMNH2": {
        "FMNH2 [reduced Flavin mononucleotide]": -1,
        "Hydrogen [H(+)]": 2,
        "FMN [oxidized Flavin mononucleotide]": 1,
    },
    "Q8H2": {"Ubiquinone-8": -1, "Hydrogen [H(+)]": 2, "Ubiquinol-8": 1},
    "MQL8": {"Menaquinone-8": -1, "Hydrogen [H(+)]": 2, "Menaquinol-8": 1},
    "DMMQL8": {
        "2-Demethylmenaquinone-8": -1,
        "Hydrogen [H(+)]": 2,
        "2-Demethylmenaquinol-8": 1,
    },
    "ACCOA": {
        "Acetyl-CoA": -1,
        "Water [H2O]": -1,
        "Hydrogen [H(+)]": 1,
        "Acetate [Acetic acid]": 1,
        "Coenzyme A": 1,
    },
    "GLU": {
        "D-Glucose": -1,
        "Water [H2O]": -1,
        "2-Oxoglutarate [Oxoglutaric acid]": 1,
        "Ammonia": 1,
        "Hydrogen [H(+)]": 2,
    },
    "PROTON": {"Hydrogen [H(+)]": 1, "Hydrogen [H(+)] transported": -1},
}  #: :meta:

EGC_SCORING_MATRIX = {"MR": 1, "RB": 3, "RF": 3, "RM": 6}  #: :meta:

################################################################################
# classes
################################################################################


class EGCSolver:
    """Parent class for the EGC solvers with generally useful functions and
    attributes. Can only be used to find, not solve EGCs directly.

    Attributes:
        - theshold: Float describing the cutoff, under which the model
            will no longer considered to be growing.
            Defaults to the MIN_GROWTH_THRESHOLD set in the growth module.
        - limit: Sets the maximal number of cores to be used.
            Defaults to 2.
        - chunksize: Chunksize to use for multiprocessing.
            Defaults to 1.
    """

    def __init__(
        self,
        threshold: float = MIN_GROWTH_THRESHOLD,
        limit: int = 2,
        chunksize: int = 1,
    ) -> None:

        # 'biological' attributes
        self.threshold = threshold

        # purely computational attributes
        self.limit = limit
        self.chunksize = chunksize

    # general functions
    # -----------------
    def check_metab_integration(
        self,
        metabolites: dict[str:int],
        model: cobra.Model,
        metab_info: Medium,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> Union[None, dict]:
        """Check if the metabolites of a reactions are in the model.
        If yes, return the dictionary of metabolites (their IDs in the model) to the factors.
        If no, return None

        Args:
            - metabolites (dict[str: int]):
                Metabolites mapped to factors.
            - model (cobra.Model):
                The model loaded with COBRApy.
            - metab_info (Medium):
                Information about the metabolites from the database,
                in for of a Medium object.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            (1) Case: metabolites not found
                    None:
                        nothing to return

            (2) Case: found
                    dict:
                        The mapping of IDs of the metabolites to the factors.
        """

        found_ids = {}

        c1_metab = [_.id for _ in model.metabolites if _.compartment == compartment[0]]
        c2_metab = [_.id for _ in model.metabolites if _.compartment == compartment[1]]

        for meta in list(metabolites.keys()):
            # get metabolite database annotations
            pos_ids = metab_info.substance_table[
                metab_info.substance_table["name"] == meta
            ][["db_type", "db_id"]]
            # check namespace availability
            if not any(namespace in _ for _ in list(pos_ids["db_type"])):
                return None
            # check metabolite availability in model
            id = pos_ids[pos_ids["db_type"].str.contains(namespace)]["db_id"]
            for i in id:
                match namespace:
                    # BiGG ID needs the compartment with a '_' as suffix
                    case "BiGG":
                        if len(metabolites) == 2:  # special case proton
                            i_p = i + "_" + compartment[1]
                            i_c = i + "_" + compartment[0]
                            if i_p not in c2_metab or i_c not in c1_metab:
                                return None
                            else:
                                found_ids[i_p] = metabolites[
                                    "Hydrogen [H(+)] transported"
                                ]
                                found_ids[i_c] = metabolites["Hydrogen [H(+)]"]
                                return found_ids
                        else:
                            i += "_" + compartment[0]
                    # No namespacce given
                    case _:
                        mes = f"Unknown namespace or no namespace given: {namespace}"
                        warnings.warn(mes)
                        return None

                if i not in c1_metab:
                    return None
                else:
                    found_ids[i] = metabolites[meta]
                    break

        return found_ids

    def add_DISSIPATIONRXNS(
        self,
        model: cobra.Model,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> cobra.Model:
        """Add the dissipation reactions a model.

        Args:
            - model (cobra.Model):
                A model loaded with COBRApy.
            - namespace (Literal['BiGG'], optional): Namespace of the model.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            cobra.Model:
                The edited model.
        """

        # retrieve information about dissipation reaction metabolites
        metab_info = Medium("dissipation reaction metabolites")
        metab_info.add_subset("DiReM")

        # add dissipations reactions to model
        for name, metabolites in DISSIPATION_RXNS.items():
            if "DISSI_" + name not in model.reactions:
                ids_in_namespace = self.check_metab_integration(
                    metabolites,
                    model,
                    metab_info,
                    namespace=namespace,
                    compartment=compartment,
                )
                if ids_in_namespace:
                    rea = cobra.Reaction(name=name, id="DISSI_" + name)
                    model.add_reactions([rea])
                    rea.add_metabolites(ids_in_namespace)

        return model

    def limit_bounds(self, model: cobra.Model):
        """Limits upper and lower bounds of

        - exchange reactions to (0, 0)
        - reversible reactions to (-1, 1)
        - irreversible reactions to (0, 1)

        Excludes dissipation reactions.

        Args:
            - model (cobra.Model):
                COBRApy model
        """
        external_comp = cobra.medium.find_external_compartment(model)
        # set fluxes for each reaction within model
        for rea in model.reactions:
            # except dissipation reactions
            if "DISSI_" in rea.id:
                continue
            # turn off exchange reactions
            elif cobra.medium.is_boundary_type(rea, "exchange", external_comp):
                rea.bounds = (0.0, 0.0)
            # limit reversible reactions to [-1, 1] -> flux 0.1
            elif rea.reversibility:
                rea.bounds = (-1.0, 1.0)
            # limit irreversible reactions to [0, 1] -> flux 0.1
            else:
                rea.bounds = (0.0, 1.0)

    # find EGCs
    # ---------

    def find_egcs(
        self,
        model: cobra.Model,
        with_reacs: bool = False,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> Union[list, tuple]:
        """Find the EGCs in a model - if exsistend.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy.
            - with_reacs (bool, optional):
                Option to either only return the names
                of the found EGC or additionally also the reactions, which
                show fluxes during testing.
                Defaults to False.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            (1) Case: ``with_reacs = False``
                    list:
                        List of found EGC names.

            (2) Case: ``with_reacs = True``
                    tuple:
                        tuple of (1) dictionary & (2) list:

                        (1) dict: dictionary of the EGCs
                        (2) list: their reactions that showed fluxes and the objective values of the test.
        """
        # add 15 energy dissipation reactions
        with model as mod_model:

            # ensure no model modifications
            mod_model = self.add_DISSIPATIONRXNS(mod_model, namespace, compartment)

            # set fluxes for each reaction within model
            self.limit_bounds(mod_model)

            # set each dissipation reaction as objective & optimize
            # options: just return the names of the found cycles
            if not with_reacs:
                egcs = []
                for name in DISSIPATION_RXNS.keys():
                    rea_id = "DISSI_" + name
                    if rea_id in mod_model.reactions:
                        mod_model.objective = "DISSI_" + name
                        solution = mod_model.optimize()

                        if (
                            solution.objective_value > MIN_GROWTH_THRESHOLD
                        ):  # optimization > 0 --> EGC detected
                            egcs.append(name)

                return egcs
            # option 2: report cycles + problematic reactions with fluxes
            else:
                egc_reactions = {}
                obj_vals = {}
                for name in DISSIPATION_RXNS.keys():
                    rea_id = "DISSI_" + name
                    if rea_id in mod_model.reactions:
                        mod_model.objective = "DISSI_" + name
                        solution = mod_model.optimize()
                        fluxes = solution.fluxes
                        objval = solution.objective_value

                        if (
                            objval > MIN_GROWTH_THRESHOLD
                        ):  # optimization > 0 --> EGC detected
                            obj_vals[name] = objval
                            egc_reactions[name] = {}

                            for rea, flux in fluxes.items():
                                # cutoff for flux
                                if (
                                    flux > MIN_GROWTH_THRESHOLD
                                    or flux < -1.0 * MIN_GROWTH_THRESHOLD
                                ) and not rea.startswith("DISSI"):
                                    egc_reactions[name][rea] = flux

                return (obj_vals, egc_reactions)

    # Helper functions for testing during solving process
    # ---------------------------------------------------

    def egcs_removed(
        self,
        model: cobra.Model,
        starting_egcs: dict,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> list:
        """Compare a list of previously found EGCs to the current
        EGCs in the model.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy after
                a try of solving the EGCs.
            - starting_egcs (dict):
                List of EGCs before trying to solve them.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            list:
                List of newly removed EGCs.
        """

        current_egcs = self.find_egcs(
            model, namespace=namespace, compartment=compartment
        )

        removed = [_ for _ in starting_egcs.keys() if _ not in current_egcs]
        added = [_ for _ in current_egcs if _ not in starting_egcs.keys()]

        if len(added) > 0:

            return []
        else:
            return removed


class GreedyEGCSolver(EGCSolver):
    """EGC solver that finds a good solution (greedy) based on modifications
    to single reactions.

    Workflow:

    - identify existing EGCs
    - test, if EGCs can be solved using single modifications of reactions

        - possible modifications:

            - deletion (RM)
            - set reversible (MR)
            - remove backward (forward only) (RB)
            - remove forward (backward only) (RF)

    - find a good - not optimal - combination of reactions, that solve
      the maximum number of EGCs that can be solved this way
    - apply solution to the model
    - report remaining EGCs, score and reactions used for solution

    Attributes:
        - all attributes of the base class :py:class:`refinegems.classes.egcs.EGCSolver`
        - scoring_matrix:
            Dictionary of the changes (RM, MR, RF, RB) against
            Integers describing the penalty scores.
    """

    def __init__(self, scoring_matrix: dict = EGC_SCORING_MATRIX, **kwargs) -> None:
        super().__init__(**kwargs)
        self.scoring_matrix = scoring_matrix

    # find modifications to solve single EGCs
    # ---------------------------------------

    def check_egc_growth(
        self,
        reac: cobra.Reaction,
        model: cobra.Model,
        bounds: tuple,
        starting_egcs: dict,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> Union[list, None]:
        """Check EGC removal and growth of a model when chaning the bounds
        of a single reaction.

        Args:
            - reac (cobra.Reaction):
                The reaction to change
            - model (cobra.Model):
                The model (COBRApy) to manipulate.
            - bounds (tuple):
                The new reactions bounds.
            - starting_egcs (dict):
                Dict of the original EGCs found in the model.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional): List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            (1) Case if EGCs removed
                    list:
                        List of EGCs that can be removed with the change

            (2) Case no removal possible
                    None:
                        no return
        """

        # set new reaction bounds
        reac.bounds = bounds
        # check egcs removal and growth
        egc_test = self.egcs_removed(model, starting_egcs, namespace, compartment)
        res = model.optimize()
        if (
            len(egc_test) > 0
            and res.status == "optimal"
            and res.objective_value > self.threshold
        ):
            return egc_test
        else:
            None

    def test_modifications(
        self,
        reaction: cobra.Reaction,
        model: cobra.Model,
        present_egc: dict,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> dict:
        """Tries four cases for a Reaction

        1. if reaction is not reversible -> make reaction reversible (MR)
        2. limit backward reaction (RB)
        3. limit forward reaction (RF)
        4. "delete" reaction by setting fluxes to 0 (RM)

        -> for each case the EGCs which are present in the model are checked if they are removed

        -> if EGCs are removed we check if the model still grows on optimal medium

        => When both limitations are True reaction is saved to corresponding dictionary


        Args:
            - reaction (cobra.Reaction):
                Reaction from a cobra.Model
            - model (cobra.Model):
                The corresponding GEM loaded with cobrapy
            - present_egc (dict):
                Dictionary of present EGCs {"egc": {}} -> EGCs are keys
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            dict:
                {"egc": {"MR":[potential_solutions],
                "RB":[potential_solutions],
                "RF":[potential_solutions],
                "RM":[potential_solutions]}}
        """
        results = {
            "ATP": {"MR": [], "RB": [], "RF": [], "RM": []},
            "CTP": {"MR": [], "RB": [], "RF": [], "RM": []},
            "GTP": {"MR": [], "RB": [], "RF": [], "RM": []},
            "UTP": {"MR": [], "RB": [], "RF": [], "RM": []},
            "ITP": {"MR": [], "RB": [], "RF": [], "RM": []},
            "NADH": {"MR": [], "RB": [], "RF": [], "RM": []},
            "NADPH": {"MR": [], "RB": [], "RF": [], "RM": []},
            "FADH2": {"MR": [], "RB": [], "RF": [], "RM": []},
            "FMNH2": {"MR": [], "RB": [], "RF": [], "RM": []},
            "Q8H2": {"MR": [], "RB": [], "RF": [], "RM": []},
            "MQL8": {"MR": [], "RB": [], "RF": [], "RM": []},
            "DMMQL8": {"MR": [], "RB": [], "RF": [], "RM": []},
            "ACCOA": {"MR": [], "RB": [], "RF": [], "RM": []},
            "GLU": {"MR": [], "RB": [], "RF": [], "RM": []},
            "PROTON": {"MR": [], "RB": [], "RF": [], "RM": []},
        }

        # skip exchange reactions
        if reaction.boundary:
            return None

        # test modifications of the reaction
        with model as m:
            original_bounds = reaction.bounds
            # case 1: make reaction reversible
            if not reaction.reversibility:  # skip if reaction is already reversible
                res = self.check_egc_growth(
                    reaction,
                    m,
                    (-1000.0, 1000.0),  # irreversible -> reversible
                    present_egc,
                    namespace,
                    compartment,
                )
                if res:
                    for egc in res:
                        results[egc]["MR"].append(reaction.id)  # MR = make reversible

            # case 2: limit backward reaction
            if original_bounds != (0.0, 1000.0):
                res = self.check_egc_growth(
                    reaction,
                    m,
                    (0.0, 1000.0),  # Remove Backward (RB) reaction
                    present_egc,
                    namespace,
                    compartment,
                )
                if res:
                    for egc in res:
                        results[egc]["RB"].append(reaction.id)

            # case 3: limit forward reaction
            if original_bounds != (-1000.0, 0):
                res = self.check_egc_growth(
                    reaction,
                    m,
                    (-1000.0, 0),  # Remove Forward (RF) reaction
                    present_egc,
                    namespace,
                    compartment,
                )
                if res:
                    for egc in res:
                        results[egc]["RF"].append(reaction.id)

            # case 4: "delete" reaction
            if original_bounds != (0.0, 0.0):
                res = self.check_egc_growth(
                    reaction,
                    m,
                    (0.0, 0.0),  # Remove (RM) reaction
                    present_egc,
                    namespace,
                    compartment,
                )
                if res:
                    for egc in res:
                        results[egc]["RM"].append(reaction.id)

        return results

    def find_mods_resolve_egcs_greedy(
        self,
        model: cobra.Model,
        present_egcs: dict,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> dict:
        """Find the (single) modifications to reactions in a cobra.Model and returns these in a dictionary.
        Splits the modification check in multiple processes.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy.
            - present_egcs (dict):
                Dict of the original EGCs found in the model.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            dict:
                Dictionary of potential modifications to resolve EGCs
                {"egc": {"MR":[potential_solutions],
                "RB":[potential_solutions],
                "RF":[potential_solutions],
                "RM":[potential_solutions]}}
        """

        output_list = []

        if len(present_egcs) == 0:
            print("No EGCs present, nothing to solve.")

        else:
            print("_____________________________________")
            print(
                f"Try to resolve the following EGCs: {[egc for egc in present_egcs]} \nThis might take a while..."
            )

            # try to find BOF
            pos_bofs = test_biomass_presence(model)
            if pos_bofs:
                bof = pos_bofs[0]
            else:
                mes = "No growth or biomass objectuve function in model. Cannot solve EGCs."
                raise KeyError(mes)

            with model as m:
                set_bounds_to_default(m)
                m.objective = bof

                # might limit processes with Pool(process=limit) -> otherwise it consumes all cores
                # limit to half of machine cores -> at least for me no speed increment with more cores
                try:
                    pool = Pool(processes=self.limit)

                    # partial -> creates function with fixed variables to call with each iteration
                    # needed since pool.imap cannot do that
                    part_test_mods = partial(
                        self.test_modifications,
                        model=m,
                        present_egc=present_egcs,
                        namespace=namespace,
                        compartment=compartment,
                    )
                    # increment in chunksize will reduce computation time -> but progressbar update also...
                    for res in list(
                        tqdm(
                            pool.imap_unordered(
                                func=part_test_mods,
                                iterable=m.reactions,
                                chunksize=self.chunksize,
                            ),
                            total=len(m.reactions),
                            desc="Resolve EGCs",
                        )
                    ):
                        if res:
                            output_list.append(res)
                finally:
                    pool.close()
                    pool.join()

            # merge output_list to the final results
            results = {}
            for output_dict in output_list:
                for egc, modifications in output_dict.items():
                    for mod, reactions in modifications.items():
                        if not egc in results.keys():
                            results[egc] = {}
                        if mod in results[egc].keys():
                            results[egc][mod] = list(set(results[egc][mod] + reactions))
                        else:
                            results[egc][mod] = list(set(reactions))

            return results

    # find a solution to solve all EGCs that
    # have a solution based on previous step
    # --------------------------------------
    def find_solution_greedy(self, results: dict, egc_reactions: dict) -> tuple:
        """Based on the originally found EGCs and the output of
        :py:func:`find_mods_resolve_egcs_greedy`, find a solution that
        solves all EGCs that can be solved with the results.

        Args:
            - results (dict):
                Output of :py:func:`find_mods_resolve_egcs_greedy`.
            - egc_reactions (dict):
                Output of :py:meth:`~refinegems.classes.egcs.EGCSolver.find_egcs` with 'with_reac=True'.
                Should be the EGCs before calculating any solutions and applying them.
            - scoring_matrix (dict, optional):
                Dictionary of the modifications types
                (RM, MR, RF, RB) and their penality score.
                Defaults to the in-build scoring matrix.

        Returns:
            tuple:
                Tuple of (1) dict & (2) int:

                (1) dictionary of reaction IDs and their mode of change
                (2) score of the solution
        """

        solved_egcs = set()
        score = 0
        reacs_for_solu = {}

        # make table and drop EGCs that are not in the model
        solution_table = pd.DataFrame(results)[list(egc_reactions.keys())]
        for col in solution_table.columns:
            solution_table[col] = solution_table[col].apply(
                lambda x: ", ".join(map(str, x))
            )
        # find not solvable (with the current code) cyles
        not_solvable = []
        for col in egc_reactions.keys():
            if all(len(_) == 0 for _ in solution_table[col]):
                not_solvable.append(col)
        # report the problems for manual curation
        if len(not_solvable) > 0:
            print(
                f"The following EGCs cannot be fixed with the current code: {not_solvable}"
            )
            solution_table.drop(not_solvable, axis=1, inplace=True)
        # fix the EGCs that can be fixed with single changes
        solution_table = solution_table.T
        for egc in solution_table.index.to_list():
            # check if sgc has already been taken care of
            if egc in solved_egcs:
                continue

            # best solution for this cycle
            for mode in self.scoring_matrix.keys():
                if solution_table.loc[egc, mode]:
                    reac = solution_table.loc[egc, mode].split(", ")[0]
                    newly_solved = [
                        _
                        for _ in solution_table.index[
                            solution_table[mode].str.contains(reac)
                        ].to_list()
                        if _ not in solved_egcs
                    ]
                    if len(newly_solved) > 0:
                        solved_egcs = set([*solved_egcs, *newly_solved])
                        score += self.scoring_matrix[mode]
                        reacs_for_solu[reac] = mode

        return (reacs_for_solu, score)

    # apply the found solution the model
    # ----------------------------------
    def apply_modifications(self, model: cobra.Model, solution: dict):
        """Apply the modifications to reactions in solution to the model.

        4 modifications are possible:

        - "RM" -> removes the reaction
        - "RB" -> removes the backwards reaction
        - "RF" -> removes the forward reaction
        - "MR" -> makes reaction reversible

        Args:
            - model (cobra.Model):
                Input model
            - solution (dict):
                Best solution from calculation in `py:func:find_solution_greedy`.
        """
        for reac, mode in solution.items():
            reaction = model.reactions.get_by_id(reac)
            if type(reaction) is cobra.Reaction:
                match mode:
                    case "RM":
                        reaction.delete()
                    case "RB":
                        reaction.bounds = (0.0, 1000.0)
                    case "RF":
                        model = reaction.bounds = (1000.0, 0.0)
                    case "MR":
                        model = reaction.bounds = (1000.0, 1000.0)
                    case _:
                        mes = (
                            f"{mode} is no viable modification... Something went wrong."
                        )
                        warnings.warn(mes)

    # run the complete solving process
    # --------------------------------
    def solve_egcs(
        self,
        model: cobra.Model,
        namespace: Literal["BiGG"] = "BiGG",
        compartment: list = ["c", "e"],
    ) -> Union[dict, None]:
        """Run the complete greedy EGC solving process.

        Note: The input model gets changed, if EGCs can be solved.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy.
            - namespace (Literal['BiGG'], optional):
                String for the namespace used in the model.
                Current options include 'BiGG'.
                Defaults to 'BiGG'.
            - compartment (list, optional):
                List of length 2 with the names of the
                compartments for the dissipations reactions.
                Defaults to ['c','e'].

        Returns:
            dict:
                Dictionary with the following entries.

                - 'solution': List of reactions for the solution.
                - 'score': Score of the solution.
                - 'remaining egcs': List of EGCs that could not be solved.

        """
        egc_reactions, obj_vals = self.find_egcs(
            model, with_reacs=True, namespace=namespace, compartment=compartment
        )
        results = self.find_mods_resolve_egcs_greedy(
            model=model,
            present_egcs=egc_reactions,
            namespace=namespace,
            compartment=compartment,
        )
        if results:
            solution, score = self.find_solution_greedy(results, egc_reactions)
            self.apply_modifications(model, solution)
            still_egc = self.find_egcs(
                model, namespace=namespace, compartment=compartment
            )

            return {"solution": solution, "score": score, "remaining egcs": still_egc}
