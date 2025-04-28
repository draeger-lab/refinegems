#!/usr/bin/env python
"""Provides functions to load and write models, media definitions and the manual annotation table

Depending on the application the model needs to be loaded with COBRApy (e.g. memote)
or with libSBML (e.g. activation of groups). Some might even require both (e.g. gap filling).
The manual_annotations table has to follow the specific layout given in the data folder in order to work with this module.
"""

__author__ = "Carolin Brune, Tobias Fehrenbach, Famke Baeuerle and Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################


import cobra
import gffutils
import logging
import os
import pandas as pd
import re
import sqlalchemy
import sqlite3

from ols_client import EBIClient
from Bio import SeqIO
from libsbml import Model as libModel
from libsbml import SBMLReader, writeSBMLToFile, SBMLValidator, SBMLDocument
from pathlib import Path
from typing import Literal, Union

from .databases import PATH_TO_DB

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# models
# ------


def load_model(
    modelpath: Union[str, list[str]], package: Literal["cobra", "libsbml"]
) -> Union[cobra.Model, list[cobra.Model], libModel, list[libModel]]:
    """Load a model.

    Args:
        - modelpath (str | list[str]):
            Path to the model or list of paths to models (string format).
        - package (Literal['cobra','libsbml']):
            Package to use to load the model.

    Returns:
        cobra.Model|list[cobra.Model]|libModel|list[libModel]:
            The loaded model(s).
    """

    def load_cobra_model(modelpath: str) -> cobra.Model:
        """Load a model using COBRApy.

        Args:
            modelpath (str): Path to the model.
                Can be a xml, json, yml or mat file.

        Raises:
            - ValueError: Unknown file extension

        Returns:
            cobra.Model:
                The loaded model object.
        """
        extension = os.path.splitext(modelpath)[1].replace(".", "")

        match extension:
            case "xml" | "sbml":
                data = cobra.io.read_sbml_model(modelpath)
            case "json":
                data = cobra.io.load_json_model(modelpath)
            case "yml" | "yaml":
                data = cobra.io.load_yaml_model(modelpath)
            case "mat":
                data = cobra.io.load_matlab_model(modelpath)
            case _:
                raise ValueError("Unknown file extension for model: ", extension)

        return data

    def load_libsbml_model(modelpath: str) -> libModel:
        """Load a model with libsbml.

        Args:
            modelpath (str): Path to the model. Should be xml.

        Returns:
            libModel:
                The loaded model object
        """

        extension = os.path.splitext(modelpath)[1].replace(".", "")

        match extension:
            case "xml" | "sbml":
                reader = SBMLReader()
                read = reader.readSBMLFromFile(modelpath)  # read from file
                mod = read.getModel()
            case "json" | "yml" | "mat":
                data = load_cobra_model(modelpath)
                mod = cobra.io.sbml._model_to_sbml(data).getModel()
            case _:
                raise ValueError("Unknown file extension for model: ", extension)

        return mod

    match modelpath:
        # read in multiple models
        case list():

            loaded_models = []
            for m in modelpath:
                if package == "cobra":
                    loaded_models.append(load_cobra_model(m))
                elif package == "libsbml":
                    loaded_models.append(load_libsbml_model(m))
                else:
                    mes = f"Unknown type for package. Unable to load model."
                    raise ValueError(mes)
            return loaded_models

        # read in a single model
        case str():

            if package == "cobra":
                return load_cobra_model(modelpath)
            elif package == "libsbml":
                return load_libsbml_model(modelpath)
            else:
                mes = f"Unknown type for package. Unable to load model."
                raise ValueError(mes)

        case Path():

            modelpath = str(modelpath)
            if package == "cobra":
                return load_cobra_model(modelpath)
            elif package == "libsbml":
                return load_libsbml_model(modelpath)
            else:
                mes = f"Unknown type for package. Unable to load model."
                raise ValueError(mes)

        case _:
            mes = f"Unkown type for modelpath: {modelpath}"
            raise TypeError(mes)


def load_document_libsbml(modelpath: str) -> SBMLDocument:
    """Loads model document using libSBML

    Args:
        - modelpath (str):
            Path to GEM

    Returns:
        SBMLDocument:
            Loaded document by libSBML
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    return read


def write_model_to_file(
    model: Union[libModel, cobra.Model], filename: Union[str, Path]
):
    """Save a model into a file.

    Args:
        - model (libModel|cobra.Model):
            The model to be saved
        - filename (str|Path):
            The filename/path to save the model to.

    Raises:
        - ValueError: Unknown file extension for model
        - TypeError: Unknown model type
    """

    def _write_cobra_model_to_file(model: cobra.Model, filepath: Path):
        """Helper function of :py:func:`~refingems.utility.io.write_model_to_file`.
        Write a model loaded with COBRApy to a file.

        Args:
            - model (cobra.Model):
                The COBRApy model.
            - filepath (Path):
                The file path to save the model to. Extensions sets the file format.

        Raises:
            - ValueError: Unknown file extension for model
        """
        try:
            extension = filepath.suffix.replace(".", "")
            match extension:
                case "xml":
                    cobra.io.write_sbml_model(model, filepath)
                case "json":
                    cobra.io.save_json_model(model, filepath)
                case "yml":
                    cobra.io.save_yaml_model(model, str(filepath))
                case "mat":
                    cobra.io.save_matlab_model(model, filepath)
                case _:
                    raise ValueError("Unknown file extension for model: ", extension)
            logging.info("Modified model written to " + str(filepath))
        except OSError as e:
            print("Could not write to file. Wrong path?")

    # Cast filename to Path object if string is provided
    if isinstance(filename, str):
        filename = Path(filename)

    # save cobra model
    if isinstance(model, cobra.core.model.Model):
        _write_cobra_model_to_file(model, filename)

    # save libsbml model
    elif isinstance(model, libModel):
        try:
            extension = filename.suffix.replace(".", "")
            match extension:
                case "xml":
                    new_document = model.getSBMLDocument()
                    writeSBMLToFile(new_document, str(filename))
                case "json" | "yml" | "mat":
                    data = cobra.io.sbml._sbml_to_model(model)
                    _write_cobra_model_to_file(data, filename)
                case _:
                    raise ValueError("Unknown file extension for model: ", extension)
            logging.info("Modified model written to " + str(filename))
        except OSError as e:
            print("Could not write to file. Wrong path?")
    # unknown model type or no model
    else:
        message = f"Unknown model type {type(model)}. Cannot save."
        raise TypeError(message)


# media
# -----


def load_substance_table_from_db(
    mediumname: str, database: str, type: Literal["testing", "standard"] = "standard"
) -> pd.DataFrame:
    """Load a substance table from a database.

    Currently available types:

    - 'testing': for debugging
    - 'standard': The standard format containing all information in long format.

    Note: 'documentation' currently object to change

    Args:
        - name (str):
            The name (or identifier) of the medium.
        - database (str):
            Path to the database.
        - type (Literal['testing','standard'], optional):
            How to load the table. Defaults to 'standard'.

    Raises:
        - ValueError: Unknown type for loading the substance table.

    Returns:
        pd.DataFrame:
            The substance table in the specified type retrieved from the database.
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    match type:
        # use this for debugging
        case "testing":
            result = cursor.execute(
                """SELECT substance.id, substance.name, substance.formula, medium2substance.flux 
                                    FROM medium, medium2substance, substance
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id
                                    """,
                (mediumname,),
            )
            substance_table = result.fetchall()
            substance_table = pd.DataFrame(
                substance_table, columns=["id", "name", "formula", "flux"]
            )

        # create table with all information, standard for generating the Medium table
        case "standard":
            result = cursor.execute(
                """SELECT substance.name, substance.formula, medium2substance.flux , medium2substance.source, substance2db.db_id, substance2db.db_type
                                    FROM medium, medium2substance, substance, substance2db
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id AND substance2db.substance_id = substance.id
                                    """,
                (mediumname,),
            )
            substance_table = result.fetchall()
            substance_table = pd.DataFrame(
                substance_table,
                columns=["name", "formula", "flux", "source", "db_id", "db_type"],
            )

        # default: throw error
        case _:
            connection.close()
            raise ValueError(f"Unknown type for loading the substance table: {type}")

    # close connection
    connection.close()

    return substance_table


def load_subset_from_db(subset_name: str) -> tuple[str, str, pd.DataFrame]:
    """Load a subset from the database.

    Args:
        - subset_name(str):
            Name of the subset to be loaded.

    Returns:
        tuple of (1) str, (2) str and (3) pd.DataFrame
            (1) name of the subset
            (2) description of the subset
            (3) substance table for the subset
    """

    # open connection to database
    connection = sqlite3.connect(PATH_TO_DB)
    cursor = connection.cursor()

    # check if subset name is valid
    result = cursor.execute(
        """SELECT * FROM subset WHERE subset.name = ?""", (subset_name,)
    )
    check = result.fetchone()
    if check:
        # set name
        name = subset_name
        # set description
        description = check[2]
        # retrieve subset from database
        db_res = cursor.execute(
            """SELECT substance.name, subset2substance.percent
                                FROM substance, subset, subset2substance
                                WHERE subset.name = ? AND subset.id = subset2substance.subset_id AND subset2substance.substance_id = substance.id
                                """,
            (subset_name,),
        )
        substance_table = pd.DataFrame(db_res.fetchall(), columns=["name", "percent"])

    return (name, description, substance_table)


# other
# -----


def load_a_table_from_database(
    table_name_or_query: str, query: bool = True
) -> pd.DataFrame:
    """| Loads the table for which the name is provided or a table containing all rows for which the query evaluates to
       | true from the refineGEMs database ('data/database/data.db')

    Args:
        - table_name_or_query (str):
            Name of a table contained in the database 'data.db'/ a SQL query
        - query (bool):
            Specifies if a query or a table name is provided with table_name_or_query

    Returns:
        pd.DataFrame:
            Containing the table for which the name was provided from the database 'data.db'
    """
    table_name_or_query = (
        sqlalchemy.text(table_name_or_query) if query else table_name_or_query
    )
    sqlalchemy_engine_input = f"sqlite:///{PATH_TO_DB}"
    engine = sqlalchemy.create_engine(sqlalchemy_engine_input)
    open_con = engine.connect()

    db_table = pd.read_sql(table_name_or_query, open_con)

    open_con.close()
    return db_table


def parse_dict_to_dataframe(str2list: dict) -> pd.DataFrame:
    """| Parses dictionary of form {str: list} &
       | Transforms it into a table with a column containing the strings and a column containing the lists

    Args:
        str2list (dict):
            Dictionary mapping strings to lists

    Returns:
        pd.DataFrame:
            Table with column containing the strings and column containing the lists
    """
    # Get max number of list length
    max_len_of_list = max(map(len, str2list.values()))

    # Fill lists with None until all lists have the same size -> Required for pd.DataFrame
    for key in str2list:
        current_list = str2list.get(key)
        while len(current_list) != max_len_of_list:
            str2list.get(key).append(None)

    df = pd.DataFrame.from_dict(str2list).stack().T.reset_index()
    df = df.drop("level_0", axis=1)

    return df


def validate_libsbml_model(model: libModel) -> int:
    """Debug method: Validates a libSBML model with the libSBML validator

    Args:
        - model (libModel):
            A libSBML model

    Returns:
        int:
            Integer specifying if validate was successful or not
    """
    validator = SBMLValidator()
    doc = model.getSBMLDocument()

    return validator.validate(doc)


# FASTA
# -----
def create_missing_genes_protein_fasta(
    fasta: str,
    missing_genes: pd.DataFrame,
    outdir: str = None,
) -> str:
    """Creates a FASTA file containing proteins for missing_genes

    .. note::

        Please keep in mind that the input FASTA file has to have Genbank format.

    Args:
        - fasta (str):
            Path to the FASTA protein file.
        - missing_genes (pd.DataFrame):
            The table of missing genes.
        - outdir (str, optional):
            Path to a directory to write the output to.
            Defaults to None.

    Returns:
        str:
            Path to the FASTA protein file for the missing genes.
    """

    # format the missing genes' locus tags
    locus_tags_dict = {
        _: "[locus_tag=" + _ + "]" for _ in list(missing_genes["locus_tag"])
    }

    # parse the protein FASTA
    protfasta = SeqIO.parse(fasta, "fasta")

    # extract the sequences of the missing genes only
    missing_seqs = []
    for seq in protfasta:
        # Case 1: locus tag equals FASTA header
        if seq.id in locus_tags_dict.keys():
            missing_seqs.append(seq)
        # Case 2: locus tag as descriptor
        elif any([k for k, v in locus_tags_dict.items() if v in seq.description]):
            for k, v in locus_tags_dict.items():
                if v in seq.description:
                    seq.id = k
                    missing_seqs.append(seq)
                    break
        # Case _: locus tag either in model or errounous
        else:
            pass

    # save the collected sequences in a new file
    if outdir:
        outfile = str(Path(outdir, "missing_genes.fasta"))
    else:
        outfile = str(Path("missing_genes.fasta"))
    SeqIO.write(missing_seqs, outfile, "fasta")

    return outfile


# GFF
# ---

def parse_gff_for_cds(
    gffpath: str, keep_attributes: dict[str:str] = None
    ) -> Union[pd.DataFrame, tuple[pd.DataFrame, str]]:
    """Parses a GFF file to obtain a mapping for the corresponding attributes listed in keep_attributes

    Args:
        - gffpath (str):
            Path to the GFF file.
        - keep_attributes (dict, optional):
            Dictionary of attributes to be kept and the corresponding column for the table.
            Defaults to None.

    Returns:
        (1) Case: ``return_variety = False``
                pd.DataFrame:
                    Dataframe containing a mapping for the corresponding attributes listed in keep_attributes


        (2) Case: ``return_variety = True``
                tuple:
                    tuple of (1) pd.DataFrame & (2) str:

                    (1) pd.DataFrame: Dataframe containing a mapping for the corresponding attributes listed in keep_attributes
                    (2) str: Found variety for provided GFF
    """

    # load the gff
    gff = gffutils.create_db(gffpath, ":memory:", merge_strategy="create_unique")
    # extract the attributes of the CDS
    cds = pd.DataFrame.from_dict([_.attributes for _ in gff.features_of_type("CDS")])
    if "locus_tag" in cds.columns:
        cds = cds.explode("locus_tag")
    genes = pd.DataFrame.from_dict([_.attributes for _ in gff.features_of_type("gene")])
    if "locus_tag" in genes.columns:
        genes = genes.explode("locus_tag")
    if "old_locus_tag" in genes.columns and "locus_tag" in genes.columns:
        cds = cds.merge(
            genes[["locus_tag", "old_locus_tag"]], how="left", on="locus_tag"
        )
        cds.drop(columns=["locus_tag"], inplace=True)
        cds.rename(columns={"old_locus_tag": "locus_tag"}, inplace=True)
    # keep only certain columns
    if keep_attributes:
        cds = cds[[_ for _ in cds.columns if _ in list(keep_attributes.keys())]]
        # rename columns
        cds.rename(
            columns={k: v for k, v in keep_attributes.items() if k in cds.columns},
            inplace=True,
        )
    if "locus_tag" in cds.columns:
        cds = cds.explode("locus_tag")

    return cds


# GBFF
# ----


def parse_gbff_for_cds(file_path: str) -> pd.DataFrame:
    """Retrieves a table containg information about the following qualifiers from a
    Genbank file: ['protein_id','locus_tag','db_xref','old_locus_tag','EC_number'].

    Args:
        - file_path (str):
            Path to the Genbank (.gbff) file.

    Returns:
        pd.DataFrame:
            A table containing the information above.
            Has the following  columns= ['ncbi_accession_version', 'locus_tag_ref','old_locus_tag','GeneID','EC number'].
    """

    temp_table = pd.DataFrame(
        columns=[
            "ncbi_accession_version",
            "locus_tag_ref",
            "old_locus_tag",
            "GeneID",
            "EC number",
        ]
    )
    attributes = ["protein_id", "locus_tag", "old_locus_tag", "db_xref", "EC_number"]

    for record in SeqIO.parse(file_path, "genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    temp_list = []
                    for a in attributes:
                        if a in feature.qualifiers.keys():
                            temp_list.append(feature.qualifiers[a][0])
                        else:
                            temp_list.append("-")
                    temp_table.loc[len(temp_table)] = temp_list

    # reformat
    pat = re.compile(r"\D")
    temp_table["GeneID"] = [pat.sub("", x) for x in temp_table["GeneID"]]

    return temp_table


# else:
# -----


def search_sbo_label(sbo_number: str) -> str:
    """Looks up the SBO label corresponding to a given SBO Term number

    Args:
        - sbo_number (str):
            Last three digits of SBO-Term as str

    Returns:
        str:
            Denoted label for given SBO Term
    """
    sbo_number = str(sbo_number)
    client = EBIClient()
    sbo = client.get_term("sbo", "http://biomodels.net/SBO/SBO_0000" + sbo_number)
    return sbo["_embedded"]["terms"][0]["label"]
