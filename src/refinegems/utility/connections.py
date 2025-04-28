#!/usr/bin/env python
"""Provides functions / connections to other tools for easier access and usage."""

__author__ = "Famke Baeuerle und Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import json
import memote
import pandas as pd
import shutil
import subprocess
import tempfile
import time
import warnings

from BOFdat import step1
from BOFdat import step2
from BOFdat.util import update
from BOFdat.util.update import determine_coefficients

from importlib.resources import files
from libsbml import Model as libModel
from MCC import MassChargeCuration
from pathlib import Path
from libsbml import readSBML
from sboannotator.SBOannotator import sbo_annotator
from typing import Literal, Union

from memote.support import consistency

# needed by memote.support.consistency
from memote.support import consistency_helpers as con_helpers

from .util import test_biomass_presence, is_stoichiometric_factor
from .io import write_model_to_file

# note:
#    for BOFdat to run correctly, one needs to change 'solution.f' to 'solution.objective_value'
#    in the coenzymes_and_ions.py file of BOFdat
#    -> see forked version of it in the Draeger-lab github

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# BOFdat
# ------


def adjust_BOF(
    genome: str,
    model_file: str,
    model: cobra.Model,
    dna_weight_fraction: float,
    weight_frac: float,
) -> str:
    """Adjust the model's BOF using BOFdat. Currently implemented are step 1
    DNA coefficients and step 2.

    Args:
        - genome (str):
            Path to the genome (e.g. .fna) FASTA file.
        - model_file (str):
            Path to the sbml (.xml) file of the model.
        - model (cobra.Model):
            The genome-scale metabolic model (from the string above), loaded with COBRApy.
        - dna_weight_fraction (float):
            DNA weight fraction for BOF step 1.
        - weight_frac (float):
            Weight fraction for the second step of BOFdat (enzymes and ions)

    Returns:
        str:
            The updated BOF reaction as a reaction string.
    """

    # BOFdat step 1:
    # --------------
    # dna coefficients
    dna_coefficients = step1.generate_dna_coefficients(
        genome, model_file, DNA_WEIGHT_FRACTION=dna_weight_fraction
    )
    bd_step1 = {}
    for m in dna_coefficients:
        bd_step1[m.id] = dna_coefficients[m]

    # BOFdat step 2:
    # --------------
    # find inorganic ions
    selected_metabolites = step2.find_coenzymes_and_ions(model_file)
    # determine coefficients
    bd_step2 = determine_coefficients(selected_metabolites, model, weight_frac)
    bd_step2.update(bd_step1)

    # update BOF
    # ----------
    #  retrieve previously used BOF
    growth_func_list = test_biomass_presence(model)
    if len(growth_func_list) == 1:
        objective_list = model.reactions.get_by_id(growth_func_list[0]).reaction.split(
            " "
        )
    elif len(growth_func_list) > 1:
        mes = f"Multiple BOFs found. Using {growth_func_list[0]} for BOF adjustment."
        warnings.warn(mes, category=UserWarning)
        objective_list = model.reactions.get_by_id(growth_func_list[0]).reaction.split(
            " "
        )
    # else not needed, as if there is no BOF, a new one will be created

    # ...............................................................
    objective_reactant = {}
    objective_product = {}
    product = False
    # get reactants, product and factors from equation
    for s in objective_list:
        if s == "+":
            continue
        elif ">" in s:
            product = True
        elif is_stoichiometric_factor(s):
            factor = s
        else:
            if product:
                objective_product[s] = factor
            else:
                objective_reactant[s] = factor

    # update BOF information with data from BOFdat
    for m in bd_step2:
        if bd_step2[m] < 0:
            objective_reactant[m] = bd_step2[m] * (-1)
        else:
            objective_product[m] = bd_step2[m]

    # create new objective function
    new_objective = " + ".join(
        "{} {}".format(value, key) for key, value in objective_reactant.items()
    )
    new_objective += " --> "
    new_objective += " + ".join(
        "{} {}".format(value, key) for key, value in objective_product.items()
    )

    return new_objective


# DIAMOND
# -------


def run_DIAMOND_blastp(
    fasta: str,
    db: str,
    sensitivity: Literal[
        "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"
    ] = "more-sensitive",
    coverage: float = 95.0,
    threads: int = 2,
    outdir: str = None,
    outname: str = "DIAMOND_blastp_res.tsv",
) -> str:
    """Run DIAMOND in BLASTp mode.

    Args:
        - fasta (str):
            The FASTA file to BLAST for.
        - db (str):
            The DIAMOND database file to BLAST against
        - sensitivity (Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], optional):
            Sensitivity mode for DIAMOND.
            Defaults to 'more-sensitive'.
        - coverage (float, optional):
            A parameter for DIAMOND
            Coverage theshold for the hits.
            Defaults to 95.0.
        - threads (int, optional):
            A parameter for DIAMOND.
            Number of threds to be used while BLASTing.
            Defaults to 2.
        - outdir (str, optional):
            Path to a directory to write the output files to.
            Defaults to None.
        - outname (str, optional):
            Name of the result file (name only, not a path).
            Defaults to 'DIAMOND_blastp_res.tsv'.

    Returns:
        str:
            Path to the results of the DIAMOND BLASTp run.
    """

    if outdir:
        outname = str(Path(outdir, "DIAMOND_blastp_res.tsv"))
        logfile = Path(outdir, "log_DIAMOND_blastp.txt")
    else:
        outname = str(Path(outname))
        logfile = Path("log_DIAMOND_blastp.txt")

    completed_blast = subprocess.run(
        [
            "diamond",
            "blastp",
            "-d",
            db,
            "-q",
            fasta,
            "--" + sensitivity,
            "--query-cover",
            str(coverage),
            "-p",
            str(threads),
            "-o",
            outname,
            "--outfmt",
            str(6),
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
        shell=False,
        stderr=subprocess.PIPE,
        text=True,
    )
    with open(logfile, "a") as f:
        f.write(completed_blast.stderr)

    return outname


def filter_DIAMOND_blastp_results(
    blasttsv: str, pid_theshold: float = 90.0
) -> pd.DataFrame:
    """Filter the results of a DIAMOND BLASTp run (see
    :py:func:`~refinegems.utility.connections.run_DIAMOND_blastp`)
    by percentage identity value (PID) and extract the matching pairs of query
    and subject IDs.

    Args:
        - blasttsv (str):
            Path to the DIAMOND BLASTp result file.
        - pid_theshold (float, optional):
            Threshold value for the PID. Given in percent.
            Defaults to 90.0.

    Raises:
        - ValueError: PID threshold has to be between 0.0 and 100.0

    Returns:
        pd.DataFrame:
            A table with the columns query_ID and subject_ID containing hits from
            BLAST run with s PID higher than the given threshold value.
    """

    if pid_theshold > 100.0 or pid_theshold < 0.0:
        raise ValueError("PID threshold has to be between 0.0 and 100.0")

    # load diamond results
    diamond_results = pd.read_csv(blasttsv, sep="\t", header=None)
    diamond_results.columns = [
        "query_ID",
        "subject_ID",
        "PID",
        "align_len",
        "no_mismatch",
        "no_gapopen",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "E-value",
        "bitscore",
    ]
    # filter by PID
    diamond_results = diamond_results[diamond_results["PID"] >= pid_theshold]
    # trim cols
    diamond_results = diamond_results[["query_ID", "subject_ID"]]

    return diamond_results


# MCC - MassChargeCuration
# ------------------------


def perform_mcc(model: cobra.Model, dir: str, apply: bool = True) -> cobra.Model:
    """Run the MassChargeCuration toll on the model and optionally directly apply
    the solution.

    Args:
        - model (cobra.Model):
            The model to use the tool on.
        - dir (str):
            Path of the directory to save MCC output in.
        - apply (bool, optional):
            If True, model is directly updated with the results.
            Defaults to True.

    Returns:
        cobra.Model:
            The model (updated or not)
    """

    # make temporary directory to save files for MCC in
    with tempfile.TemporaryDirectory() as temp:

        # use MCC
        if apply:
            # update model
            balancer = MassChargeCuration(model, update_ids=False, data_path=temp)
        else:
            # do not change original model
            model_copy = model.copy()
            balancer = MassChargeCuration(model_copy, update_ids=False, data_path=temp)

    # save reports
    balancer.generate_reaction_report(Path(dir, model.id + "_mcc_reactions"))
    balancer.generate_metabolite_report(Path(dir, model.id + "_mcc_metabolites"))
    balancer.generate_visual_report(Path(dir, model.id + "_mcc_visual"))

    return model


# Memote
# -------


def run_memote(
    model: cobra.Model,
    type: Literal["json", "html"] = "html",
    return_res: bool = False,
    save_res: Union[str, None] = None,
    verbose: bool = False,
) -> Union[dict, str, None]:
    """Run the memote snapshot function on a given model loaded with COBRApy.

    Args:
        - model (cobra.Model):
            The model loaded with COBRApy.
        - type (Literal['json','html'], optional):
            Type of report to produce.
            Can be 'html' or 'json'.
            Defaults to 'html'.
        - return_res (bool, optional):
            Option to return the result.
            Defaults to False.
        - save_res (str | None, optional):
            If given a path string, saves the report
            under the given path. Defaults to None.
        - verbose (bool, optional):
            Produce a more verbose ouput.
            Defaults to False.

    Raises:
        - ValueError: Unknown input for parameter type

    Returns:
        (1) Case ``return_res = True`` and ``type = json``:
                dict:
                    The json dictionary.

        (2) Case ``return_res = True`` and ``type = html``:
                str:
                    The html string.

        (3) Case ``return_res = False``:
                None:
                    no return
    """

    # verbose output I
    if verbose:
        print("\n# -------------------\n# Analyse with MEMOTE\n# -------------------")
        start = time.time()

    # run memote
    ret, res = memote.suite.api.test_model(
        model,
        sbml_version=None,
        results=True,
        pytest_args=None,
        exclusive=None,
        skip=None,
        experimental=None,
        solver_timeout=10,
    )

    # load depending on type
    match type:
        case "html":
            snap = memote.suite.api.snapshot_report(res, html=True)
            result = snap
        case "json":
            snap = memote.suite.api.snapshot_report(res, html=False)
            result = json.loads(snap)
        case _:
            message = f"Unknown input for parameter type: {type} "
            raise ValueError(message)

    # option to save report
    if save_res:
        with open(save_res, "w", encoding="utf-8") as f:
            f.write(result)

    # verbose output II
    if verbose:
        end = time.time()
        print(f"\ttotal time: {end - start}s")

    # option to return report
    if return_res:
        return result


def get_memote_score(memote_report: dict) -> float:
    """Extracts MEMOTE score from report

    Args:
        - memote_report (dict):
            Output from :py:func:`~refinegems.utility.connections.run_memote`.

    Returns:
        float:
            MEMOTE score
    """
    return memote_report["score"]["total_score"]


# SBOannotator
# ------------


def run_SBOannotator(model: libModel) -> libModel:
    """Run SBOannotator on a model to annotate the SBO terms.

    Args:
        - model (libModel):
            The model loaded with libSBML.

    Returns:
        libModel:
            The model with corrected / added SBO terms.
    """

    dbs_scheme = files("sboannotator").joinpath("create_dbs.sql")

    with tempfile.TemporaryDirectory() as tempdir:
        write_model_to_file(model, str(Path(tempdir, "tempmodel.xml")))
        # run SBOannotator
        doc = readSBML(str(Path(tempdir, "tempmodel.xml")))
        model = doc.getModel()
        copy_scheme = shutil.copy(dbs_scheme, Path(tempdir, "dbs.sql"))
        model = sbo_annotator(
            doc,
            model,
            "constrained-based",
            str(Path(tempdir, "dbs")),
            str(Path(tempdir, "dud.xml")),
        )
    return model
