#!/usr/bin/env python
"""Access information from different databases or compare a model or model entities
with them. This module provides variables and function for accessing databases
for better model curation and annotation.

The following databases have functionalities implemented:

- BiGG
- ChEBI
- KEGG
- ModelSEED
- NCBI
- UniProt
"""

__author__ = """Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune, 
            Jan-Philipp Leusch and Tobias Fehrenbach"""

############################################################################
# requirements
############################################################################

import cobra
import io
import libchebipy
import numpy as np
import pandas as pd
import re
import requests
import warnings
import xmltodict

from Bio import Entrez, SeqIO
from Bio.KEGG import REST, Gene, Enzyme
from bioservices.kegg import KEGG
from cobra import Model as cobraModel
from typing import Literal, Union
from tqdm import tqdm

tqdm.pandas()
pd.options.mode.chained_assignment = None  # suppresses the pandas SettingWithCopyWarning; comment out before developing!!

from .connections import run_DIAMOND_blastp, filter_DIAMOND_blastp_results
from .io import load_a_table_from_database, create_missing_genes_protein_fasta

############################################################################
# variables
############################################################################

# database urls
# -------------
BIGG_METABOLITES_URL = "http://bigg.ucsd.edu/api/v2/universal/metabolites/"  #: :meta:

# Compartments in BiGG namespace
# ------------------------------
ALL_BIGG_COMPARTMENTS_ONE_LETTER = (
    "c",
    "e",
    "p",
    "m",
    "x",
    "r",
    "v",
    "n",
    "g",
    "u",
    "l",
    "h",
    "f",
    "s",
    "i",
    "w",
    "y",
)  #: :meta:
ALL_BIGG_COMPARTMENTS_TWO_LETTER = ("im", "cx", "um", "cm", "mm")  #: :meta:

# BioCyc databases prefixes
# -------------------------
BIOCYC_TIER1_DATABASES_PREFIXES = ["META", "ECO", "ECOLI", "HUMAN"]  #: :meta:

############################################################################
# functions
############################################################################

# BiGG
# ----
"""Map an ID to information in BiGG or compare model entities to BiGG.
"""


def get_BiGG_metabs_annot_via_dbid(
    metabolite: cobra.Metabolite, id: str, dbcol: str, compartment: str = "c"
) -> None:
    """Search for a BiGG ID and add it to a metabolite annotation.
    The search is based on a column name of the BiGG metabolite table
    and an ID to search for. Additionally, using the given compartment name,
    the found IDs are filtered for matching compartments.

    Args:
        - metabolite (cobra.Metabolite):
            The metabolite. Needs to a a COBRApy Metabolte object.
        - id (str):
            The ID to search for in the database.
        - dbcol (str):
            Name of the column of the database to check the ID against.
        - compartment (str, optional):
            The compartment name. Needs to be a valid BiGG compartment ID.
            Defaults to 'c'.
    """
    # check, if a BiGG annotation is given
    if not "bigg.metabolite" in metabolite.annotation.keys():
        # seach the database for the found id
        bigg_search = load_a_table_from_database(
            f"SELECT * FROM bigg_metabolites WHERE '{dbcol}' = '{id}'", query=True
        )
        if len(bigg_search) > 0:
            # check, if the matches also match the compartment
            metabolite.annotation["bigg.metabolite"] = [
                _.rsplit("_", 1)[0]
                for _ in bigg_search["id"].tolist()
                if _.endswith(f"_{compartment}")
            ]
            # add final matches to the annotations of the metabolite
            if len(metabolite.annotation["bigg.metabolite"]) == 0:
                metabolite.annotation.pop("bigg.metabolite")


def add_annotations_from_BiGG_metabs(metabolite: cobra.Metabolite) -> None:
    """Check a cobra.metabolite for bigg.metabolite annotations. If they exists,
    search for more annotations in the BiGG database and add them to the metabolite.

    Args:
        - metabolite (cobra.Metabolite):
            The metabolite object.
    """
    if "bigg.metabolite" in metabolite.annotation.keys():
        bigg_information = load_a_table_from_database(
            "SELECT * FROM bigg_metabolites WHERE id = '"
            + f"' OR id = '".join(metabolite.annotation["bigg.metabolite"])
            + "'",
            query=True,
        )
        db_id_bigg = {
            "BioCyc": "biocyc",
            "MetaNetX (MNX) Chemical": "metanetx.chemical",
            "SEED Compound": "seed.compound",
            "CHEBI": "chebi",
            "KEGG Compound": "kegg.compound",
            "InChI Key": "inchikey",
        }
        for db in db_id_bigg:
            info = list(set(bigg_information[db].dropna().to_list()))
            if len(info) > 0:
                info = ",".join(info)
                info = [
                    x.strip() for x in info.split(",")
                ]  # make sure all entries are a separate list object
                if db_id_bigg[db] in metabolite.annotation.keys():
                    metabolite.annotation[db_id_bigg[db]] = list(
                        set(info + metabolite.annotation[db_id_bigg[db]])
                    )
                else:
                    metabolite.annotation[db_id_bigg[db]] = info


def _add_annotations_from_bigg_reac_row(row: pd.Series, reac: cobra.Reaction) -> None:
    """Given a row of the BiGG reaction database table and a cobra.Reaction object,
    extend the annotation of the latter with the information of the former.

    Args:
        - row (pd.Series):
            The row of the database table.
        - reac (cobra.Reaction):
            The reaction object.
    """

    dbnames = {
        "RHEA": "rhea",
        "BioCyc": "biocyc",
        "MetaNetX (MNX) Equation": "metanetx.reaction",
        "EC Number": "ec-code",
        "SEED Reaction": "seed.reaction",
        "Reactome Reaction": "reactome",
        "KEGG Reaction": "kegg.reaction",
    }
    for dbname, dbprefix in dbnames.items():
        if row[dbname]:
            ids_to_add = row[dbname].strip().split(",")
            if dbprefix in reac.annotation.keys():
                reac.annotation[dbprefix] = list(
                    set(reac.annotation[dbprefix]).union(set(ids_to_add))
                )
            else:
                reac.annotation[dbprefix] = ids_to_add


# ChEBI
# -----
"""Primarily used for curating the media database. Obtain information about 
metabolites (table format)."""


def add_info_from_ChEBI_BiGG(
    missing_metabs: pd.DataFrame, charge=True, formula=True, iupac=True
) -> pd.DataFrame:
    """Adds information from CHEBI/BiGG to the provided dataframe.

    The following informations can be added:

    - charge
    - formula
    - iupac (name)

    Args:
       - missing_metabs (pd.DataFrame):
          Table containing metabolites & the respective ChEBI & BiGG IDs

    Returns:
       pd.DataFrame:
          Input table extended with the charges & chemical formulas obtained from ChEBI/BiGG.
    """

    # check if a row contains a ChEBI ID, take the first and make sure its in the format: CHEBI:234567
    def get_chebi_id(row: pd.Series) -> str:

        chebi = row.get("ChEBI")
        if pd.isnull(chebi):
            return None
        elif (
            type(chebi) == str
        ):  # check for mulitple entries (assuming first one is the best fitting one)
            chebi = chebi.split(",")[0]
            if "CHEBI" in chebi:
                return chebi
            else:
                return "CHEBI:" + chebi
        else:
            return "CHEBI:" + str(chebi)

    # Finds the charges through the ChEBI/BiGG API, defaults to: 0
    def find_charge(row: pd.Series) -> int:
        chebi_id = get_chebi_id(row)
        bigg_id = str(row.get("bigg_id"))
        charge = None
        if chebi_id:  # Get charge from ChEBI (Returns always a charge)
            chebi_entity = libchebipy.ChebiEntity(chebi_id)
            return chebi_entity.get_charge()
        elif bigg_id != "nan":  # Get charge from BiGG if no ChEBI ID available
            try:
                charge = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()[
                    "charges"
                ][
                    0
                ]  # Take first charge
            except ValueError:
                pass
            # If no charge was found, charge=0
            return charge if charge else 0

    # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
    def find_formula(row: pd.Series) -> str:
        chebi_id = get_chebi_id(row)
        bigg_id, chem_form = str(row.get("bigg_id")), str(row.get("Chemical Formula"))
        chem_formula = None
        if chebi_id:  # Get formula from ChEBI
            chebi_entity = libchebipy.ChebiEntity(chebi_id)
            chem_formula = chebi_entity.get_formula()
        if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
            if bigg_id != "nan":  # Get formula from BiGG
                try:
                    chem_formula = requests.get(
                        BIGG_METABOLITES_URL + bigg_id[:-2]
                    ).json()["formulae"][
                        0
                    ]  # Take first formula
                except ValueError:
                    pass
            if not chem_formula:  # If no formula was found with BiGG ID
                # Get formula already existing in dataframe or set to 'No formula'
                chem_formula = chem_form if chem_form != "nan" else "No formula"
        return chem_formula

    # using the ChEBI ID, retrieve the IUPAC name from the ChEBI database.
    def find_iupac(row: pd.Series) -> str:

        chebi_id = get_chebi_id(row)
        if chebi_id:
            chebi_entity = libchebipy.ChebiEntity(chebi_id)
            # only take the IUPAC names
            iupac_names = []
            for name in chebi_entity.get_names():
                print(name)
                if name.get_source() == "IUPAC":
                    iupac_names.append(name.get_name())
                    break

            if len(iupac_names) > 0:
                iupac_name = ", ".join(iupac_names)
            else:
                iupac_name = None
        else:
            iupac_name = None

        return iupac_name

    if charge:
        missing_metabs["charge"] = missing_metabs.apply(find_charge, axis=1)
    if formula:
        missing_metabs["New Chemical Formula"] = missing_metabs.apply(
            find_formula, axis=1
        )
        missing_metabs["Chemical Formula"] = missing_metabs["New Chemical Formula"]
        missing_metabs.drop("New Chemical Formula", axis=1, inplace=True)
    if iupac:
        missing_metabs["ChEBI_IUPAC"] = missing_metabs.apply(find_iupac, axis=1)

    return missing_metabs


# KEGG
# ----
"""Retrieve and parse information from the KEGG database. 
"""


def get_kegg_genes(organismid: str) -> pd.DataFrame:
    """Extracts list of genes from KEGG given an organism

    Args:
        - organismid (str):
            KEGG ID of organism which the model is based on

    Returns:
        pd.DataFrame:
            Table of all genes denoted in KEGG for the organism
    """
    k = KEGG()
    gene_list = k.list(organismid)

    return pd.read_table(io.StringIO(gene_list), header=None)


def parse_KEGG_gene(locus_tag: str) -> dict:
    """Based on a locus tag, fetch the corresponding KEGG entry and
    parse it into a dictionary containing the following information (if available):

    - ec-code
    - orthology
    - references

    Args:
        - locus_tag (str):
            The locus in the format <organism_id>:<locus_tag>

    Returns:
        dict:
            The collected information.
    """

    gene_info = dict()
    gene_info["orgid:locus"] = locus_tag

    # retireve KEGG gene entry
    try:
        gene_entry = list(Gene.parse(REST.kegg_get(locus_tag)))[0]
    except Exception as e:
        gene_entry = None

    # skip, if no entry found
    if not gene_entry:
        gene_info["ec-code"] = None
        return gene_info

    # extract orthology and ec-code
    if len(gene_entry.orthology) > 0:
        # get KEGG orthology ID
        kegg_orthology = [_[0] for _ in gene_entry.orthology]
        gene_info["kegg.orthology"] = kegg_orthology
        # get EC number
        ec_numbers = [
            re.search(r"(?<=EC:).*(?=\])", orth[1]).group(0)
            for orth in gene_entry.orthology
            if re.search(r"(?<=EC:).*(?=\])", orth[1])
        ]
        if isinstance(ec_numbers, list) and len(ec_numbers) > 0:
            gene_info["ec-code"] = [
                ec for ec_str in ec_numbers for ec in ec_str.split(" ")
            ]

    if not "ec-code" in gene_info.keys():
        gene_info["ec-code"] = None

    # get more information about connections to other databases
    if len(gene_entry.dblinks) > 0:
        for dbname, ids in gene_entry.dblinks:
            conform_dbname = re.sub(
                pattern=r"(NCBI)(.*)(ID$)", repl="\\1\\2", string=dbname
            )  # Remove ID if NCBI in name
            conform_dbname = re.sub(
                r"[^\w]", "", conform_dbname
            )  # remove special signs except underscore
            conform_dbname = conform_dbname.lower()  # make lower case
            gene_info[conform_dbname] = ids

    return gene_info


def parse_KEGG_ec(ec: str) -> dict:
    """Based on an EC number, fetch the corresponding KEGG entry and
    parse it into a dictionary containing the following information (if available):

    - ec-code
    - id (kegg.reference)
    - equation
    - reference
    - pathway

    Args:
        - ec (str):
            The EC number in the format 'x.x.x.x'

    Returns:
        dict:
            The collected information about the KEGG entry.
    """

    ec_info = dict()
    ec_info["ec-code"] = ec

    # retrieve KEGG entry
    try:
        ec_entry = list(Enzyme.parse(REST.kegg_get(ec)))[0]
    except Exception as e:
        ec_entry = None
        ec_info["id"] = None
        ec_info["equation"] = None
        ec_info["reference"] = None
        return ec_info

    # retrieve reaction information from entry
    rn_numbers = [
        re.search(r"(?<=RN:).*(?=\])", reac).group(0)
        for reac in ec_entry.reaction
        if re.search(r"(?<=RN:).*(?=\])", reac)
    ]
    if len(rn_numbers) > 0:
        ec_info["id"] = [_.split(" ") for _ in rn_numbers]
        ec_info["equation"] = None
    else:
        ec_info["id"] = None
        ec_info["equation"] = ec_entry.reaction

    # retrieve database links from entry
    refs = dict()
    # orthology not possible with biopython
    if len(ec_entry.pathway) > 0:
        refs["kegg.pathway"] = [_[1] for _ in ec_entry.pathway]
    if len(ec_entry.dblinks) > 0:
        for dbname, ids in ec_entry.dblinks:
            if "BRENDA" in dbname:
                refs["brenda"] = ids
            if "CAS" == dbname:
                refs["cas"] = ids
    ec_info["reference"] = refs

    return ec_info


def kegg_reaction_parser(rn_id: str) -> dict:
    """Get the entry of a KEGG reaction ID and
    parse the information into a dictionary.

    Args:
        - rn_id (str):
            A reaction ID existing in KEGG.

    Returns:
        dict:
            The KEGG entry information as a dictionary.
    """

    # get KEGG reaction entry
    try:
        kegg_reac = REST.kegg_get(f"rn:{rn_id}")
        kegg_reac = kegg_reac.read()
    except Exception as e:
        return None

    # parse the entry for necessary information
    features = {"db": []}
    collect = False
    for line in kegg_reac.split("\n"):
        if line:
            if line.startswith("NAME"):
                features["name"] = line.replace("NAME", "", 1).strip()
            elif line.startswith("EQUATION"):
                features["equation"] = line.replace("EQUATION", "", 1).strip()
            elif line.startswith("ENZYME"):
                features["ec-code"] = re.split(
                    r"\s+", line.replace("ENZYME", "", 1).strip()
                )
                collect = "ec-code"
            elif line.startswith("RCLASS"):
                features["kegg.rclass"] = [
                    line.replace("RCLASS", "", 1).strip().split(" ")[0]
                ]
                collect = "kegg.rclass"
            elif line.startswith("PATHWAY"):
                features["kegg.pathway"] = [
                    line.replace("PATHWAY", "", 1).strip().split(" ")[0]
                ]
                collect = "kegg.pathway"
            elif line.startswith("DBLINKS"):
                features["db"] = [line.replace("DBLINKS", "", 1).strip()]
                collect = "db"
            elif collect and line[0] != "/":
                if len(features["db"]) == 0:
                    if line[0].isupper():
                        collect = None
                    else:
                        line = line.strip()
                        match collect:
                            case "ec-code":
                                features[collect].append(line)
                            case _:
                                features[collect].append(line.split(" ")[0])
                else:
                    features["db"].append(line.strip())
            else:
                continue

    # parse references
    db_entries = features["db"]
    features["db"] = {"kegg.reaction": rn_id}
    for entry in db_entries:
        db, identifier = entry.split(":")
        db = db.strip().lower()
        if db in features["db"]:
            features["db"][db] = features["db"][db].append(identifier)
        else:
            features["db"][db] = [identifier.strip()]

    if "ec-code" in features:
        features["db"]["ec-code"] = features["ec-code"]
        del features["ec-code"]
    if "kegg.rclass" in features:
        features["db"]["kegg.rclass"] = features["kegg.rclass"]
        del features["kegg.rclass"]
    if "kegg.pathway" in features:
        features["db"]["kegg.pathway"] = features["kegg.pathway"]
        del features["kegg.pathway"]

    return features


# ModelSEED
# ---------
""" Reports mismatches in charges and formulae based on ModelSEED

Extracts ModelSEED data from a given tsv file, extracts all metabolites from a given model. Both lists of metabolites are compared by charge and formula.
"""


def get_modelseed_compounds() -> pd.DataFrame:
    """Extracts compounds from ModelSEED which have BiGG Ids

    Returns:
        pd.DataFrame:
            Table containing ModelSEED data
    """
    # Get only rows where BiGG is contained
    com = load_a_table_from_database(
        "SELECT *, INSTR(aliases, 'BiGG:') bigg FROM modelseed_compounds WHERE bigg > 0"
    )

    def get_bigg_ids(aliases):
        try:
            aliases_list = aliases.split("|")
            bigg = [x[6:] for x in aliases_list if re.search(r"BiGG: .*", x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com["BiGG"] = com.apply(lambda row: get_bigg_ids(row["aliases"]), axis=1)

    return com.loc[:, ["id", "name", "formula", "mass", "charge", "BiGG"]]


def get_model_charges(model: cobraModel) -> pd.DataFrame:
    """Extracts all metabolites from model

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        pd.DataFrame:
            Table containing charges and formulae of model metabolites
    """
    charges = {}
    for metab in model.metabolites:
        charges[metab.id[:-2]] = [metab.charge, metab.formula]

    df_charges = (
        pd.DataFrame.from_dict(
            charges, orient="index", columns=["charge_model", "formula_model"]
        )
        .reset_index()
        .rename(columns={"index": "BiGG"})
    )

    return df_charges


def get_modelseed_charges(modelseed_compounds: pd.DataFrame) -> pd.DataFrame:
    """Extract table with BiGG, charges and formulae

    Args:
        - modelseed_compounds (pd.DataFrame):
            ModelSEED data. Output from :py:func:`~refinegems.utility.db_access.get_modelseed_compounds`.

    Returns:
        pd.DataFrame:
            Table containing charges and formulae of ModelSEED metabolites
    """
    modelseed_compounds = modelseed_compounds.loc[
        :, ["charge", "BiGG", "formula"]
    ].rename(columns={"charge": "charge_modelseed", "formula": "formula_modelseed"})
    lst_col = "BiGG"
    x = modelseed_compounds.assign(
        **{lst_col: modelseed_compounds[lst_col].str.split(";")}
    )
    df_ms = pd.DataFrame(
        {
            col: np.repeat(x[col].values, x[lst_col].str.len())
            for col in x.columns.difference([lst_col])
        }
    ).assign(**{lst_col: np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return df_ms


def compare_model_modelseed(
    model_charges: pd.DataFrame, modelseed_charges: pd.DataFrame
) -> pd.DataFrame:
    """Compares tables with charges / formulae from model & modelseed

    Args:
        - model_charges (pd.DataFrame):
            Charges and formulae of model metabolites. Output of :py:func:`~refinegems.utility.db_access.get_model_charges`.
        - modelseed_charges (pd.DataFrame):
            Charges and formulae of ModelSEED metabolites. Output of :py:func:`~refinegems.utility.db_access.get_modelseed_charges`.

    Returns:
        pd.DataFrame:
            Table containing info whether charges / formulae match
    """
    df_comp = pd.merge(model_charges, modelseed_charges, on="BiGG", how="left")

    def f(x):
        return True if float(x["charge_model"]) == x["charge_modelseed"] else False

    def g(x):
        return True if x["formula_model"] == x["formula_modelseed"] else False

    df_comp["charge_match"] = df_comp.apply(f, axis=1)
    df_comp["formula_match"] = df_comp.apply(g, axis=1)

    return df_comp


def get_charge_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with charge mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame):
            Charge and formula mismatches. Output from :py:func:`~refinegems.utility.db_access.compare_model_modelseed`.

    Returns:
        pd.DataFrame:
            Table containing metabolites with charge mismatch
    """
    return df_comp.loc[~df_comp["charge_match"]].dropna(subset=["charge_modelseed"])


def get_formula_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with formula mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame):
            Charge and formula mismatches. Output from :py:func:`~refinegems.utility.db_access.compare_model_modelseed`.

    Returns:
        pd.DataFrame:
            Table containing metabolites with formula mismatch
    """
    return df_comp.loc[~df_comp["formula_match"]].dropna(subset=["formula_modelseed"])


def get_compared_formulae(formula_mismatch: pd.DataFrame) -> pd.DataFrame:
    """Compare formula by atom pattern

    Args:
        formula_mismatch (pd.DataFrame):
            Table with column containing atom comparison. Output from :py:func:`~refinegems.utility.db_access.get_formula_mismatch`.

    Returns:
        pd.DataFrame:
            table containing metabolites with formula mismatch
    """

    def formula_comparison(f1, f2):
        # from Jan Leusch
        formula_pattern = r"[A-Z][a-z]?\\d*"
        atom_pattern = r"[A-Z][a-z]?"
        atom_number_pattern = r"\\d+"
        difference = {}
        f1_dict = {}
        f2_dict = {}
        f1 = re.findall(formula_pattern, f1)
        f2 = re.findall(formula_pattern, f2)
        for p in f1:
            key = re.findall(atom_pattern, p)[0]
            value = re.findall(atom_number_pattern, p)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f1_dict[key] = value
        for q in f2:
            key = re.findall(atom_pattern, q)[0]
            value = re.findall(atom_number_pattern, q)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f2_dict[key] = value
        difference = f1_dict
        if not f2_dict:
            return difference
        else:
            for q in f2_dict:
                if q in difference:
                    difference[q] -= f2_dict[q]
                else:
                    difference[q] = -f2_dict[q]
        for a in list(difference):
            if difference[a] == 0:
                difference.pop(a)
        return difference

    formula_mismatch["formula_comparison"] = formula_mismatch.apply(
        lambda row: formula_comparison(row["formula_model"], row["formula_modelseed"]),
        axis=1,
    )

    return formula_mismatch


def compare_to_modelseed(model: cobraModel) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Executes all steps to compare model metabolites to ModelSEED metabolites

    Args:
        - model (cobraModel):
            Model loaded with COBRApy

    Returns:
        tuple:
            Tables with charge (1) & formula (2) mismatches

            (1) pd.DataFrame: Table with charge mismatches
            (2) pd.DataFrame: Table with formula mismatches
    """
    ms_comp = get_modelseed_compounds()
    model_charges = get_model_charges(model)
    modelseed_charges = get_modelseed_charges(ms_comp)
    df_comp = compare_model_modelseed(model_charges, modelseed_charges)
    charge_mismatch = get_charge_mismatch(df_comp)
    formula_mismatch = get_formula_mismatch(df_comp)
    formula_comp = get_compared_formulae(formula_mismatch)
    return charge_mismatch, formula_comp


# NCBI
# ----
"""Connect to NBCI (using Biopython/Entrez) to extract information
like EC number of an NBCI protein accession number or information about locus tags.
"""


def _search_ncbi_for_gp(
    row: pd.Series, id_type: Literal["refseq", "ncbiprotein"]
) -> pd.Series:
    """Fetches protein name and locus tag from NCBI

    Args:
        - row (pd.Series):
            Row of a pandas DataFrame containing RefSeq/NCBI Protein IDs in columns
        - id_type (Literal['refseq', 'ncbiprotein']):
            ID type of IDs in provided row.
            Can be one of ['refseq', 'ncbiprotein'].

    Returns:
        pd.Series:
            Modified input row
    """
    match id_type:
        case "refseq":
            ncbi_id = row["REFSEQ"]
        case "ncbiprotein":
            ncbi_id = row["NCBI"]
        case _:
            raise TypeError(
                f"Provided type for id_type {id_type} invalid. Please use one of ['refseq', 'ncbiprotein']."
            )
            return row

    if not pd.isnull(ncbi_id):

        try:
            handle = Entrez.efetch(
                db="protein", id=ncbi_id, rettype="gbwithparts", retmode="text"
            )
            records = SeqIO.parse(handle, "gb")

            for i, record in enumerate(records):
                match id_type:
                    case "refseq":
                        row["name_refseq"] = record.description
                    case "ncbiprotein":
                        for feature in record.features:
                            if feature.type == "CDS":
                                row["name_ncbi"] = record.description
                                row["locus_tag"] = feature.qualifiers["locus_tag"][0]

        except Exception as e:
            print(f"{e} with {id_type} ID {ncbi_id}")

    return row


# fetching the EC number (if possible) from NCBI
# based on an NCBI protein ID (accession version number)
def get_ec_from_ncbi(mail: str, ncbiprot: str) -> Union[str, None]:
    """Based on a NCBI protein accession number, try and fetch the
    EC number from NCBI.

    Args:
        - mail (str):
            User's mail address for the NCBI ENtrez tool.
        - ncbiprot (str):
            The NCBI protein accession number.

    Returns:
        (1) Case: fetching successful
                str:
                    The EC number associated with the protein ID based on NCBI.

        (2) Case: fetching unsuccessful
                None:
                    Nothing to return
    """
    try:
        Entrez.email = mail
        handle = Entrez.efetch(db="protein", id=ncbiprot, rettype="gpc", retmode="xml")
        for feature_dict in xmltodict.parse(handle)["INSDSet"]["INSDSeq"][
            "INSDSeq_feature-table"
        ]["INSDFeature"]:
            if isinstance(feature_dict["INSDFeature_quals"]["INSDQualifier"], list):
                for qual in feature_dict["INSDFeature_quals"]["INSDQualifier"]:
                    if qual["INSDQualifier_name"] == "EC_number":
                        return qual["INSDQualifier_value"]
    except Exception as e:
        return None


# Uniprot
# -------
"""Currently mainly used for gapfilling, map your proteins against EC numbers 
using the SwissProt database"""


def map_dmnd_res_to_sp_ec_brenda(
    dmnd_results: pd.DataFrame, swissprot_mapping_path: str
) -> pd.DataFrame:
    """Map the results of a DIAMOND BLASTp run (filtered, see
    :py:func:`~refinegems.curation.db_access.db.filter_DIAMOND_blastp_results`)

    Args:
        - dmnd_results (pd.DataFrame):
            The results of the DIAMOND run.
        - swissprot_mapping_path (str):
            The path to the SwissProt mapping file (IDs against BRENDA and EC,
            for information on how to get them, refer to :py:func:`~refinegems.utility.set_up.download_url`)

    Returns:
        pd.DataFrame:
            The resulting mapping (no duplicates).
    """

    def _combine_EC_BRENDA(brenda: str, ec: str) -> list:
        """Helper function for combining the entries for BRENDA and EC number
        for one ID (SwissProt mapping file) into a list of EC numbers.

        Args:
            - brenda (str):
                The BRENDA column value.
            - ec (str):
                The EC number column value.

        Returns:
            list:
                The list with the EC numbers.
        """

        nums_list = []
        # get brenda ids
        if brenda and not pd.isna(brenda):
            for num in brenda.split(";"):
                nums_list.append(num.strip())

        # get ec numbers
        if ec and not pd.isna(ec):
            for num in ec.split(";"):
                nums_list.append(num.strip())

        # remove duplicates
        nums_list = list(set(nums_list))
        # filter non-complete EC numbers : X.X.X.X => 7 or more characters
        nums_list = [_ for _ in nums_list if len(_) >= 7]

        return nums_list

    # load the SwissProt mapping file
    swissprot_mapping = pd.read_csv(swissprot_mapping_path, sep="\t")
    swissprot_mapping.dropna(subset=["BRENDA", "EC number"], how="all", inplace=True)

    # extract SwissProts IDs from subject_ID
    dmnd_results.columns = ["locus_tag", "UniProt"]
    dmnd_results["UniProt"] = dmnd_results["UniProt"].apply(lambda x: x.split("|")[1])

    # match
    dmnd_results = dmnd_results.merge(
        swissprot_mapping, left_on="UniProt", right_on="Entry", how="left"
    )
    dmnd_results.drop("Entry", axis=1, inplace=True)
    dmnd_results["ec-code"] = dmnd_results.apply(
        lambda x: _combine_EC_BRENDA(x["BRENDA"], x["EC number"]), axis=1
    )
    dmnd_results.drop(["BRENDA", "EC number"], axis=1, inplace=True)
    dmnd_results = dmnd_results.explode("ec-code")
    dmnd_results.drop_duplicates(inplace=True)

    return dmnd_results


def get_ec_via_swissprot(
    fasta: str,
    db: str,
    missing_genes: pd.DataFrame,
    swissprot_mapping_file: str,
    outdir: str = None,
    sens: Literal[
        "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"
    ] = "more-sensitive",
    cov: float = 95.0,
    t: int = 2,
    pid: float = 90.0,
) -> pd.DataFrame:
    """Based on a protein FASTA and a missing genes tables, mapped them to EC numbers
    using a Swissprot DIAMOND database and a SwissProt mapping file (see :py:func:`~refinegems.utility.set_up.download_url`
    on how to download the needed files).

    Args:
        - fasta (str):
            Path to the FASTA protein file.
        - db (str):
            Path to the DIAMOND database (SwissProt).
        - missing_genes (pd.DataFrame):
            The table of missing genes.
        - swissprot_mapping_file (str):
            Path to the SwissProt mapping file.
        - outdir (str, optional):
            Path to a directory to write the output to.
            Defaults to None.
        - sens (Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], optional):
            Sensitivity mode of DIAMOND blastp. Defaults to 'more-sensitive'.
        - cov (float, optional):
            Coverage threshold for DIAMOND blastp. Defaults to 95.0.
        - t (int, optional):
            Number of threads to use for DIAMOND blastp. Defaults to 2.
        - pid (float, optional):
            Percentage identity value to use as a cutoff for the results
            of the DIAMOND blastp run. Defaults to 90.0.

    Returns:
        pd.DataFrame:
            The missing genes table extended by the mapping to an EC number,
            if successful.
    """

    # Step 1: Make a FASTA out of the missing genes
    miss_fasta = create_missing_genes_protein_fasta(fasta, missing_genes, outdir)
    # Step 2: Run DIAMOND
    #         blastp mode against SwissProt DB
    blast_path = run_DIAMOND_blastp(
        miss_fasta, db, sensitivity=sens, coverage=cov, threads=t, outdir=outdir
    )
    # Step 3: filter DIAMOND hits
    dmnd_res = filter_DIAMOND_blastp_results(blast_path, pid)
    # Step 4: map to Swissprot mapping file
    mapped_res = map_dmnd_res_to_sp_ec_brenda(dmnd_res, swissprot_mapping_file)
    # Step 5: Aggregate UniProt IDs for unique combinations of
    #         EC numbers and locus tags
    mapped_res = (
        mapped_res.groupby(["locus_tag", "ec-code"])
        .agg({"UniProt": lambda x: x.tolist()})
        .reset_index()
    )

    return mapped_res


# DIAMOND database
# ----------------
"""Handle and retrieve information from DIAMOND databases"""


def map_to_homologs(
    fasta: str,
    db: str,
    missing_genes: pd.DataFrame,
    mapping_file: str,
    outdir: str = None,
    sens: Literal[
        "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"
    ] = "more-sensitive",
    cov: float = 95.0,
    t: int = 2,
    pid: float = 90.0,
    email=None,
) -> pd.DataFrame:

    # Step 1: Make a FASTA out of the missing genes
    miss_fasta = create_missing_genes_protein_fasta(fasta, missing_genes, outdir)

    # Step 2: Run DIAMOND
    #         blastp mode against SwissProt DB
    blast_path = run_DIAMOND_blastp(
        miss_fasta, db, sensitivity=sens, coverage=cov, threads=t, outdir=outdir
    )
    # Step 3: filter DIAMOND hits
    dmnd_res = filter_DIAMOND_blastp_results(blast_path, pid)
    dmnd_res.rename(
        columns={"query_ID": "locus_tag", "subject_ID": "ncbiprotein"}, inplace=True
    )

    # Step 4: map to NCBI information
    # if a precomputet mapping exists
    if mapping_file:
        # load pre-compiles ncbi mapping
        ncbi_mapping = pd.read_csv(mapping_file, dtype="str")

        # merge with input info table:
        table = pd.merge(
            dmnd_res,
            ncbi_mapping,
            left_on="ncbiprotein",
            right_on="ncbi_accession_version",
        )
        table.replace({"-": None}, inplace=True)
        table.drop(columns=["ncbi_accession_version"], inplace=True, axis=1)
        table.rename(columns={"EC number": "ec-code"}, inplace=True)

        # aggregate to reduce redundant information
        table = (
            table.groupby(["locus_tag", "ec-code"], dropna=False)
            .agg(
                {
                    "ncbiprotein": "first",  # lambda x: [_ for _ in x if not pd.isna(_)]
                    "locus_tag_ref": "first",  # lambda x: [_ for _ in x if not pd.isna(_)]
                    "old_locus_tag": "first",  # lambda x: [_ for _ in x if not pd.isna(_)]
                    "GeneID": "first",  # lambda x: [_ for _ in x if not pd.isna(_)]
                }
            )
            .reset_index()
        )
        # locus_tag ncbiprotein locus_tag_ref old_locus_tag GeneID EC number

    # no mapping file
    else:
        if email:
            print(
                "\t\tNo precomputed mapping. Retrieving information directly from NCBI.\n\t\tThis may take a while."
            )
            table["locus_tag_ref"] = pd.Series(dtype="str")
            table["old_locus_tag"] = pd.Series(dtype="str")
            table["GeneID"] = pd.Series(dtype="str")
            # .......................
            # @DEBUG
            # if len(table) > 10:
            #     table = table.sample(10)
            #     print('Running in debugging mode')
            # .......................
            table["ec-code"] = table["ncbiprotein"].progress_apply(
                lambda x: get_ec_from_ncbi(x, email), axis=1
            )
        else:
            warnings.warn(
                "No email provided for NCBI quieries. Skipping mapping to EC numbers."
            )

    return table
