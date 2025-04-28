#!/usr/bin/env python
"""Provides the medium class and additional functions to handle media.

Functionalities include (amongst others):

- loading a medium into a Medium object from a database, a file or a model
- adding a medium to a model
- adding media information to the database
- extending, change and manipulate various parts of a medium to create the desired medium
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import copy
import io
import numpy as np
import pandas as pd
import random
import re
import sqlite3
import string
import sys
import warnings
import yaml

from colorama import init as colorama_init
from colorama import Fore
from pathlib import Path
from sqlite_dump import iterdump
from typing import Literal, Union, Any

from ..utility.databases import PATH_TO_DB
from ..utility.io import load_substance_table_from_db, load_subset_from_db

############################################################################
# variables
############################################################################

ALLOWED_DATABASE_LINKS = ["BiGG", "MetaNetX", "SEED", "VMH", "ChEBI", "KEGG"]  #: :meta:
REQUIRED_SUBSTANCE_ATTRIBUTES = ["name", "formula", "flux", "source"]  #: :meta:
INTEGER_REGEX = re.compile(r"^[-+]?([1-9]\d*|0)$")  #: :meta:
FLOAT_REGEX = re.compile(r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?")  #: :meta:

############################################################################
# classes
############################################################################


class Medium:
    """Class describing a medium.

    Attributes:
        - name (str):
            The name or abbreviation of the medium.
        - substance_table (pd.DataFrame):
            A table containing information about the medium in silico compounds. Long format.
        - description (str, optional):
            Short description of the medium. Defaults to None.
        - doi (str):
            Reference(s) to the original publication of the medium. Defaults to None.
    """

    def __init__(
        self,
        name: str,
        substance_table: pd.DataFrame = pd.DataFrame(
            columns=["name", "formula", "flux", "source", "db_id", "db_type"]
        ),
        description: str = None,
        doi: str = None,
    ):
        """Initialise a Medium object.

        Args:
            - name (str):
                The name or abbreviation of the medium.
            - substance_table (pd.DataFrame, optional):
                A table containing information about the medium in silico compounds. Long format.
                Defaults to an empty table with the columns ['name','formula','flux','source','db_id','db_type'].
            - description (str, optional):
                Short description of the medium. Defaults to None.
            - doi (str, optional):
                Reference(s) to the original publication of the medium.. Defaults to None.
        """

        self.name = name
        self.description = description
        self.substance_table = substance_table
        self.doi = doi

    def add_substance_from_db(self, name: str, flux: float = 10.0):
        """Add a substance from the database to the medium.

        Args:
            - name (str):
                Name of the substance. Should be part of the database substance.name column.
            - flux (float, optional):
                Sets the flux value of the new substance. Defaults to 10.0.
        """
        # build connection to DB
        connection = sqlite3.connect(PATH_TO_DB)
        cursor = connection.cursor()
        # fetch substance
        result = cursor.execute(
            """ SELECT substance.name, substance.formula, substance2db.db_id, substance2db.db_type 
                                    FROM substance, substance2db 
                                    WHERE substance.name = ? AND substance.id = substance2db.substance_id""",
            (name,),
        )
        substance = result.fetchall()

        # check if fetch was successful
        if len(substance) == 0:
            warnings.warn(f"Could not fetch substance {name} from DB.")
            return

        # add substance to table
        substance_table = pd.DataFrame(
            substance, columns=["name", "formula", "db_id", "db_type"]
        )
        substance_table.insert(2, "flux", flux)
        substance_table.insert(3, "source", None)
        self.substance_table = pd.concat(
            [self.substance_table, substance_table], ignore_index=True
        )

    def remove_substance(self, name: str):
        """Remove a substance from the medium based on its name

        Args:
            - name (str):
                Name of the substance to remove.
        """

        self.substance_table = self.substance_table[
            self.substance_table["name"] != name
        ]

    def get_source(self, element: str) -> list[str]:
        """Get the source of a given element for the medium.

        Search for the given element (elemental symbol e.g. O), excluding pattern matches that are
        followed by other lower-case letters and returm them as a list of sources for the given element.

        Args:
            - element (str):
                The symbol of the element to search the sources of

        Returns:
            list[str]:
                The list of the names of the sources (no duplicates).
        """

        return list(
            set(
                self.substance_table[
                    self.substance_table["formula"]
                    .str.contains(element + "(?![a-z])", case=True, regex=True)
                    .fillna(False)
                ]["name"]
            )
        )

    def set_source(self, element: str, new_source: str):
        """Set the source for a given element to a specific substance by deleting all
        other sources of said element before adding the new source.

        Args:
            - element (str):
                The element to set the source for, e.g. 'O' for oxygen.
            - new_source (str):
                The new source. Should be the name of a substance in the database, otherwise no new source will be set.
        """

        # get sources for element
        current_sources = self.get_source(element)
        # remove current sources from medium
        for s in current_sources:
            self.remove_substance(s)
        # add new source to medium
        self.add_substance_from_db(new_source)

    def is_aerobic(self) -> bool:
        """Check if the medium contains O2 / dioxygen.

        Returns:
            bool:
                Results of the test, True if pure oxygen is in the medium.
        """

        test = (
            self.substance_table[["name", "formula"]]
            .loc[
                (self.substance_table["name"] == "Dioxygen [O2]")
                & (self.substance_table["formula"] == "O2")
            ]
            .any()
            .all()
        )
        return test

    def make_aerobic(self, flux: float = None):
        """If the medium is curretly anaerobic, add oxygen to the medium to make it aerobic.

        Args:
            - flux(float,optional):
                The flux value for the oxygen to be added. Defaults to None.
        """
        if not self.is_aerobic():
            # build connection to DB
            connection = sqlite3.connect(PATH_TO_DB)
            cursor = connection.cursor()
            # fetch oxygen
            result = cursor.execute(
                """ SELECT substance.name, substance.formula, substance2db.db_id, substance2db.db_type 
                                    FROM substance, substance2db 
                                    WHERE substance.formula = 'O2' AND substance.name = 'Dioxygen [O2]' AND substance.id = substance2db.substance_id"""
            )
            oxygen = result.fetchall()
            # add oxygen to substance table
            oxygen_table = pd.DataFrame(
                oxygen, columns=["name", "formula", "db_id", "db_type"]
            )
            oxygen_table.insert(2, "flux", flux)
            oxygen_table.insert(3, "source", None)
            self.substance_table = pd.concat(
                [self.substance_table, oxygen_table], ignore_index=True
            )

    def make_anaerobic(self):
        """If the medium is currently aerobic, deletes the oxygen from it to make it anaerobic."""
        if self.is_aerobic():
            # remove dioxygen // O2 from the substance table
            self.substance_table.drop(
                self.substance_table[
                    (self.substance_table["name"] == "Dioxygen [O2]")
                    & (self.substance_table["formula"] == "O2")
                ].index,
                inplace=True,
            )

    def set_default_flux(
        self, flux: float = 10.0, replace: bool = False, double_o2: bool = True
    ):
        """Set a default flux for the model.

        Args:
            - flux (float, optional):
                Default flux for the medium. Defaults to 10.0.
            - replace (bool, optional):
                Replace al fluxes with the default.
                Defaults to False.
            - double_o2 (bool, optional):
                Tag to double the flux for oxygen only. Works only with replace=True. Defaults to True.
        """
        # replace fluxes
        if replace:
            self.substance_table["flux"] = flux
            if self.is_aerobic() and double_o2:
                self.substance_table.loc[
                    self.substance_table["name"] == "Dioxygen [O2]", "flux"
                ] = (2 * flux)
        # keep already set fluxes
        else:
            self.substance_table["flux"] = self.substance_table["flux"].fillna(flux)
            if (
                self.is_aerobic()
                and double_o2
                and self.substance_table.loc[
                    self.substance_table["name"] == "Dioxygen [O2]", "flux"
                ].any()
            ):
                self.substance_table.loc[
                    self.substance_table["name"] == "Dioxygen [O2]", "flux"
                ] = (2 * flux)

    def set_oxygen_percentage(self, perc: float = 1.0):
        """Set oxygen percentage of the medium.

        Args:
            - perc (float, optional):
                Percentage of oxygen. Defaults to 1.0 (= 100%)
        """

        if self.is_aerobic():
            self.substance_table.loc[
                self.substance_table["name"] == "Dioxygen [O2]", "flux"
            ] = (
                self.substance_table.loc[
                    self.substance_table["name"] == "Dioxygen [O2]", "flux"
                ]
                * perc
            )
        else:
            warnings.warn(
                f"WARNING: no oxygen detected in medium {self.name}, cannot set oxygen percentage."
            )

    def combine(
        self,
        other: "Medium",
        how: Union[Literal["+"], None, float, tuple[float]] = "+",
        default_flux: float = 10.0,
        inplace: bool = False,
    ) -> "Medium":
        """Combine two media into a new one.

        Modes to combine media (input for param how):

        - None -> combine media, remove all flux values (= set them to None). Sets sources to None as well.
        - '+' -> Add fluxes of the same substance together.
        - float -> Calculate flux * percentage (float) for first medium and flux * 1.0-percentage (float) for second medium and add fluxes of same substance together.
        - tuple(float,float) -> Same as above, except both percentages are given.

        Args:
            - other (Medium):
                The medium to combine with.
            - how (Union[Literal['+'],None,float,tuple[float]], optional):
                How to combine the two media.
                Options listed in header. Defaults to '+'.
            - default_flux (float, optional):
                Flux to use in combine-modes (except how=None) for NaN/None values.
                Defaults to 10.0.

        Returns:
            Medium:
                The combined medium.
        """

        def combine_media_with_fluxes(
            combined: "Medium", second_medium: "Medium"
        ) -> "Medium":
            """Helper function for :py:func:`combine`. Add a substance table to a medium by
            combining fluxes of substances of the same name.

            Args:
                - combined (Medium):
                    The medium that will be the combined one.
                - second_medium (Medium):
                    The substance table to integrate into the medium.

            Returns:
                Medium:
                    The combined medium.
            """
            # get dict of name + flux of second medium
            possible_dups = (
                second_medium[["name", "flux"]].groupby("name").first().reset_index()
            )
            possible_dups = dict(zip(possible_dups.name, possible_dups.flux))
            # search for overlapping substances and add fluxes together
            added = []
            for subs, flux in possible_dups.items():
                if subs in combined.substance_table["name"].values:
                    combined.substance_table.loc[
                        combined.substance_table["name"] == subs, "flux"
                    ] += flux
                    added.append(subs)
            # add remaining substances to medium from the second one
            second_medium = second_medium[~second_medium["name"].isin(added)]
            combined.substance_table = pd.concat(
                [combined.substance_table, second_medium], ignore_index=True
            )

            return combined

        # combine description and Co
        if inplace:
            combined = self
        else:
            combined = copy.deepcopy(self)
        combined.name = self.name + "+" + other.name
        combined.description = (
            f"Combined medium constructed from {self.name} and {other.name}."
        )
        combined.doi = str(self.doi) + ", " + str(other.doi)

        # combine substance table
        match how:

            # add fluxes together
            case "+":
                # replace None values
                combined.substance_table["flux"] = combined.substance_table[
                    "flux"
                ].fillna(default_flux)
                to_add = copy.deepcopy(other.substance_table)
                to_add["flux"] = to_add["flux"].fillna(default_flux)
                # add together
                combined = combine_media_with_fluxes(combined, to_add)

            # add fluxes based on percentages
            case float() | tuple():
                if type(how) == float:
                    how = (how, 1.0 - how)
                # replace None values
                combined.substance_table["flux"] = combined.substance_table[
                    "flux"
                ].fillna(default_flux)
                second_medium = copy.deepcopy(other.substance_table)
                second_medium["flux"] = second_medium["flux"].fillna(default_flux)
                # calculate new fluxes based on percentages
                combined.substance_table["flux"] = combined.substance_table.flux.apply(
                    lambda x: x * how[0]
                )
                second_medium["flux"] = second_medium.flux.apply(lambda x: x * how[1])
                # add together
                combined = combine_media_with_fluxes(combined, second_medium)

            # combine and set all fluxes to 0
            case None:
                combined.substance_table = pd.concat(
                    [combined.substance_table, other.substance_table], ignore_index=True
                )
                # remove fluxes and sources, as they are no longer a fit
                combined.substance_table["flux"] = None
                combined.substance_table["source"] = None
                # remove duplicate rows
                combined.substance_table.drop_duplicates(
                    inplace=True, ignore_index=True
                )

            case _:
                raise ValueError(f"Unknown input for parameter how: {how}")

        return combined

    def __add__(self, other: "Medium") -> "Medium":
        return self.combine(other)

    def add_subset(
        self, subset_name: str, default_flux: float = 10.0, inplace: bool = True
    ) -> "Medium":
        """Add a subset of substances to the medium, returning a newly generated one.

        Args:
            - subset_name (str):
                The type of subset to be added. Name should be in database-substset-id.
            - default_flux (float, optional):
                Default flux value to calculate fluxes from based  on the percentages saved in the database.
                Defaults to 10.0.

        Returns:
            Medium:
                A new medium that is the combination of the set subset and the old one.
                In the case that the given subset name is not found in the database,
                the original medium is returned.
        """

        # open connection to database
        connection = sqlite3.connect(PATH_TO_DB)
        cursor = connection.cursor()

        # check if subset name is valid
        result = cursor.execute(
            """SELECT 1 FROM subset WHERE subset.name = ?""", (subset_name,)
        )
        check = result.fetchone()
        if check:

            # retrieve subset from database
            db_res = cursor.execute(
                """SELECT substance.name, substance.formula, subset2substance.percent, substance2db.db_id, substance2db.db_type
                                    FROM substance, subset, subset2substance, substance2db
                                    WHERE subset.name = ? AND subset.id = subset2substance.subset_id AND subset2substance.substance_id = substance.id AND substance.id = substance2db.substance_id
                                    """,
                (subset_name,),
            )
            subs = db_res.fetchall()

            # reformat
            # --------
            # rename
            subs = pd.DataFrame(
                subs, columns=["name", "formula", "percent", "db_id", "db_type"]
            )
            # percent -> flux
            subs["flux"] = subs["percent"].apply(lambda x: default_flux * x)
            subs.drop("percent", axis=1, inplace=True)
            # add source column
            subs["source"] = None
            # sort
            colorder = ["name", "formula", "flux", "source", "db_id", "db_type"]
            subs = subs.reindex(columns=colorder)

            # create subset
            sub_medium = Medium(subset_name, subs, description=f"subset {subset_name}")

            # combine with current
            return self.combine(sub_medium, inplace=inplace)

        else:
            warnings.warn(
                f"Could not find subset in DB, nothing added to medium: {subset_name}"
            )
            # just return the original medium
            return self

    # functions for export table
    # --------------------------

    def produce_medium_docs_table(self, folder: str = "./", max_width: int = 80) -> str:
        """Produces a rst-file containing reStructuredText for the substance table for documentation.

        Args:
            - folder (str, optional):
                Path to folder/directory to save the rst-file to. Defaults to './'.
            - max_width (int, optional):
                Maximal table width of the rst-table. Defaults to 80.
        """

        def calculate_column_widths_docs(
            header: list, max_width: int, flux_width=8
        ) -> str:
            """Helper function for :py:func:`produce_medium_docs_table`.
            Calculates the columns widths based on the header lengths and maximal widths.

            Args:
                - header (list):
                    List of column names of the table.
                - max_width (int):
                    Maximal widths of the output table.
                - flux_width (int, optional):
                    Widths of the flux column, if 'flux' is part of 'header'. Defaults to 8.

            Returns:
                str: The columns widths as a string, separated by a whitespace.
            """

            if len(header) == 2:
                partition = max_width // 3
                return f"{str(max_width-partition)} {partition}"
            else:
                partition = (max_width - flux_width) // 2
                return f"{str(max_width-flux_width-partition)} {flux_width} {partition}"

        def produce_medium_docs_table_row(row: pd.Series, file: io.TextIOWrapper):
            """Helper function for :py:func:`produce_medium_docs_table`.
            Tranforms each row of the substance table into a row of the rst-file.

            Args:
                - row (pd.Series):
                    The row of the Medium.substance_table.
                - file (io.TextIOWrapper):
                    The connection to the file to write the rows into.
            """

            list = row.to_list()
            file.write(f"  * - {list[0]}\n")
            for l in list[1:]:
                file.write(f"    - {l}\n")

        # make sure given directory path ends with '/'
        if not folder.endswith("/"):
            folder = folder + "/"

        with open(folder + f"{self.name}.rst", "w") as f:

            # slim table to columns of interest for documentation
            m_subs = self.substance_table[["name", "flux", "source"]]
            m_subs.drop_duplicates(keep="first", inplace=True)

            if all(m_subs["flux"].values == None):
                m_subs = m_subs.drop("flux", axis=1)

            header = m_subs.columns

            widths = calculate_column_widths_docs(header, max_width)

            title = f"{self.description} ({self.name})"

            # Produce header/title of HTML page
            f.write(f"{title}\n" f"{'^' * len(title)}\n\n")

            # produce descriptor
            f.write(
                f".. list-table:: {self.name} composition\n"
                f"  :name: {self.name.lower()}_comp\n"
                "  :align: center\n"
                f"  :widths: {widths}\n"
                "  :header-rows: 1\n"
                "  :class: no-scrollbar-table\n\n"
            )

            # produce header
            f.write(f"  * - {header[0]}\n")
            for l in header[1:]:
                f.write(f"    - {l}\n")

            # produce table body
            m_subs.apply(produce_medium_docs_table_row, file=f, axis=1)

            f.close()

    def export_to_file(self, type: str = "tsv", dir: str = "./", max_widths: int = 80):
        """Export medium, especially substance table.

        Args:
            - type (str, optional):
                Type of file to export to. Defaults to 'tsv'. Further choices are 'csv', 'docs', 'rst'.
            - dir (str, optional):
                Path to the directory to write the file to. Defaults to './'.
            - max_widths (int, optional):
                Maximal table width for the documentation table (). Only viable for 'rst' and 'docs'.
                Defaults to 80.

        Raises:
            - ValueError: Unknown export type if type not in ['tsv','csv','docs','rst']
        """

        match type:
            case "tsv":
                self.substance_table.to_csv(
                    Path(dir, self.name + "_min_medium" + ".tsv"), sep="\t", index=False
                )
            case "csv":
                self.substance_table.to_csv(
                    Path(dir, self.name + "_min_medium" + ".csv"), sep=";", index=False
                )
            case "docs" | "rst":
                self.produce_medium_docs_table(folder=dir, max_width=max_widths)
            case _:
                raise ValueError("Unknown export type: {type}")

    # functions for conversion
    # ------------------------
    def export_to_cobra(
        self,
        namespace: Literal["Name", "BiGG"] = "BiGG",
        default_flux: float = 10.0,
        replace: bool = False,
        double_o2: bool = True,
    ) -> dict[str, float]:
        """Export a medium to the COBRApy format for a medium.

        Args:
            - namespace (Literal['Name', 'BiGG'], optional):
                Namespace to use. Defaults to 'BiGG'.
            - default_flux (float, optional):
                Default flux to substitute missing values. Defaults to 10.0.
            - replace (bool, optional):
                Replace all values with the default flux. Defaults to False.
            - double_o2 (bool, optional):
                Double the flux of oxygen. Defaults to True.

        Raises:
            - ValueError: Unknown namespace.

        Returns:
            dict[str,float]:
                The exported medium.
        """

        self.set_default_flux(default_flux, replace=replace, double_o2=double_o2)

        match namespace:
            case "Name":

                names = self.substance_table[["name", "flux"]]
                names["EX"] = biggs["db_id"] + " exchange"
                cobra_medium = pd.Series(names.flux.values, index=names.EX).to_dict()

            case "BiGG":

                biggs = self.substance_table[
                    self.substance_table["db_type"].str.contains("BiGG")
                ][["name", "db_id", "flux"]]
                biggs["db_id_EX"] = "EX_" + biggs["db_id"] + "_e"
                cobra_medium = pd.Series(
                    biggs.flux.values, index=biggs.db_id_EX
                ).to_dict()

            case _:
                raise ValueError(f"Unknown namespace: {namespace}")

        return cobra_medium


############################################################################
# functions for loading and writing from DB
############################################################################


def load_medium_from_db(
    name: str, database: str = PATH_TO_DB, type: str = "standard"
) -> Medium:
    """Load a medium from a database.

    Args:
        - name (str):
            The name (or identifier) of the medium.
        - database (str, optional):
            Path to the database. Defaults to the in-built database.
        - type (str, optional):
            How to load the medium. Defaults to 'standard'.

    Raises:
        - ValueError: Unknown medium name.

    Returns:
        Medium:
            The medium retrieved from the database.
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # get description
    result = cursor.execute(
        "SELECT medium.description FROM medium WHERE medium.name = ?", (name,)
    )
    if result:
        description = result.fetchone()[0]  # each name should be unique
    else:
        raise ValueError(f"Unknown medium name: {name}")

    # get DOI(s)
    result = cursor.execute(
        "SELECT medium.reference FROM medium WHERE medium.name = ?", (name,)
    )
    doi = result.fetchone()[0]

    # close connection to DB
    connection.close()

    # get substance table
    substance = load_substance_table_from_db(name, database, type)

    return Medium(
        name=name, substance_table=substance, description=description, doi=doi
    )


def generate_docs_for_subset(subset_name: str, folder: str = "./", max_width: int = 80):
    """Generate documentation for a subset.

    Args:
        - subset_name (str):
            Name of the subset.
        - folder (str, optional):
            Folder to save the output to. Defaults to './'.
        - max_width (int, optional):
            Maximal table width for the documentation page. Defaults to 80.
    """

    name, description, subs = load_subset_from_db(subset_name)

    with open(Path(folder, f"{name}.rst"), "w") as f:

        partition = max_width // 3
        width = f"{str(max_width-partition)} {partition}"
        header = subs.columns

        # Produce header/title of HTML page
        f.write(f"{name}\n" f"{'^' * len(name)}\n")

        f.write(description + "\n\n")

        # produce descriptor
        f.write(
            f".. list-table:: {name} composition\n"
            f"  :name: {name.lower()}_comp\n"
            "  :align: center\n"
            f"  :widths: {width}\n"
            "  :header-rows: 1\n"
            "  :class: no-scrollbar-table\n\n"
        )

        # produce header
        f.write(f"  * - {header[0]}\n")
        for l in header[1:]:
            f.write(f"    - {l}\n")

        # produce table body
        subs.apply(Medium.produce_medium_docs_table.produce_medium_docs_table_row, file=f, axis=1)


############################################################################
# functions for reading media from extern
############################################################################


def load_media(yaml_path: str) -> tuple[list[Medium], list[str, None]]:
    """Load the information from a media configuration file.

    Args:
        - yaml_path (str):
            The path to a media configuration file in YAML-format.

    Returns:
        tuple[list[Medium],list[str,None]]:
            Tuple of two lists (1) & (2)

            (1) list: list of the loaded media and
            (2) list: list of supplement modes
    """

    media_list = []
    supplement_list = []

    with open(yaml_path, "r") as stream:

        loaded = yaml.safe_load(stream)
        params = loaded["params"] if "params" in loaded.keys() else None
        media = loaded["media"] if "media" in loaded.keys() else None

        # handle internal/in-build media compositions
        if media:

            for name, p in media.items():

                # load base medium from either
                # base
                if p and "base" in p.keys():
                    m_name = p["base"].split(" ")
                    if len(m_name) > 1:
                        new_medium = load_medium_from_db(m_name[1])
                        new_medium.substance_table["flux"] = new_medium.substance_table[
                            "flux"
                        ].apply(
                            lambda x: x * float(m_name[0]) if type(x) == float else x
                        )
                    else:
                        new_medium = load_medium_from_db(p["base"])
                    new_medium.name = name
                # extern
                elif p and "external_base" in p.keys():
                    m_name = p["external_base"].split(" ")
                    if len(m_name) > 1:
                        new_medium = load_external_medium(m_name[1])
                        new_medium.substance_table["flux"] = new_medium.substance_table[
                            "flux"
                        ].apply(
                            lambda x: x * float(m_name[0]) if type(x) == float else x
                        )
                    else:
                        new_medium = load_external_medium(p["external_base"])
                    new_medium.name = name
                # default name
                else:
                    new_medium = load_medium_from_db(name)

                # add subsets from DB
                if p and "add_subset" in p.keys():
                    # interate over all subsets to add
                    for subset, dflux in p["add_subset"].items():
                        new_medium = new_medium.add_subset(subset, default_flux=dflux)

                # add media from DB
                if p and "add_medium" in p.keys():
                    # interate over all media to add
                    for medium, perc in p["add_medium"].items():
                        new_medium = new_medium.combine(
                            load_medium_from_db(medium), how=perc
                        )

                # from external
                if p and "add_external" in p.keys():
                    for a, perc in p["add_external"].items():
                        new_medium = new_medium.combine(
                            load_external_medium(a), how=perc
                        )

                # check anaerobic / aerobic settings
                if p and "aerobic" in p.keys():
                    if p["aerobic"]:
                        new_medium.make_aerobic()
                    else:
                        new_medium.make_anaerobic()
                else:
                    if params and "aerobic" in params.keys():
                        if params["aerobic"]:
                            new_medium.make_aerobic()
                        else:
                            new_medium.make_anaerobic()

                # set default flux
                if p and "default_flux" in p.keys():
                    new_medium.set_default_flux(p["default_flux"], replace=True)
                elif params and "default_flux" in params.keys():
                    new_medium.set_default_flux(params["default_flux"], replace=True)

                # set o2_percentage
                if p and "o2_percent" in p.keys():
                    new_medium.set_oxygen_percentage(p["o2_percent"])
                elif params and "o2_percent" in params.keys():
                    new_medium.set_oxygen_percentage(params["o2_percent"])

                # add additional substances from DB
                if p and "add_substance" in p.keys():
                    for s, f in p["add_substance"].items():
                        # read in flux
                        if not f:
                            if p and "default_flux" in p.keys():
                                f = p["default_flux"]
                            elif params and "default_flux" in params.keys():
                                f = params["default_flux"]
                        elif isinstance(f, str):
                            if "%" in f:
                                f = float(f[: f.index("%")]) / 100
                                if p and "default_flux" in p.keys():
                                    f = f * p["default_flux"]
                                elif params and "default_flux" in params.keys():
                                    f = f * params["default_flux"]
                                else:
                                    f = f * 10.0
                            else:
                                warn_string = f"Could not read in flux for {s}: {f}. \nWill be using default 10.0 instead."
                                warnings.warn(warn_string)
                                f = 10.0
                        # change medium
                        if s in new_medium.substance_table["name"].tolist():
                            # change fluxes only
                            new_medium.substance_table.loc[
                                new_medium.substance_table["name"] == s, "flux"
                            ] = f
                        else:
                            # add substance to medium
                            new_medium.add_substance_from_db(s, f)

                # supplement settings
                if p and "supplement" in p.keys():
                    supplement_list.append(p["supplement"])
                elif params and "supplement" in params:
                    supplement_list.append(params["supplement"])
                else:
                    supplement_list.append(None)

                # append medium to list
                media_list.append(new_medium)

    return (media_list, supplement_list)


def read_substances_from_file(path: str) -> pd.DataFrame:
    """Read in a TSV with substance information into a table.

    Format of the TSV (with example):
    name \t formula \t flux \t source \t X \t X | ...
    water \t H20 \t 10.0 \t .....

    X: placeholder for database names (columns filled with corresponding IDs of the substances)
    X = see ALLOWED_DATABASE_LINKS

    Args:
        - path(str):
            The path to the input file.

    Returns:
        pd.DataFrame:
            The table of substance information read from the file
    """

    # read in the table
    substance_table = pd.read_csv(path, sep="\t", comment="#")

    # validate the table
    substance_cols = substance_table.columns
    for req in REQUIRED_SUBSTANCE_ATTRIBUTES:
        if req in substance_cols:
            continue
        else:
            substance_table[req] = None

    # remove unwanted database links
    all_allowed = REQUIRED_SUBSTANCE_ATTRIBUTES + ALLOWED_DATABASE_LINKS
    for_removal = [_ for _ in substance_cols if _ not in all_allowed]
    substance_table.drop(for_removal, axis=1, inplace=True)

    # convert all NaNs into None
    substance_table = substance_table.replace({np.nan: None})

    return substance_table


def load_external_medium(how: Literal["file", "console"], **kwargs) -> Medium:
    """Read in an external medium.

    Currently available options for how

    - 'console': read in medium via the console
    - 'file': read in medium from a file, requires a 'path=str' argument to be passed.

    About the format (console, file):

    The substances have to be saved in a TSV file table (format see :py:func:`read_substances_from_file`).
    Further information for the 'file' option have to added as comments in the format: `# info_name: info`.
    Information should but do not need to contain name, description and reference.

    Args:
        - how (Literal['file','console']):
            How (or from where) the medium should be read in.
            Available options are given above.

    Raises:
        - ValueError: Unknown description for how.

    Returns:
        Medium:
            The read-in medium.
    """

    match how:
        # interactive mode using the console
        case "console":

            # prompt for name
            name = input("Enter a medium name:\n")

            # prompt for description
            description = input("Enter a short description for the medium:\n")

            # promt for description
            reference = input("Enter the reference link (DOI) for the medium:\n")

            # prompt for path to substance table
            substance_path = input("Enter a path to a table of substances (TSV):\n")
            # tsv format
            substances = read_substances_from_file(substance_path)

            # construct medium
            new_medium = Medium(
                name=name,
                description=description,
                substance_table=substances,
                doi=reference,
            )

            return new_medium

        # read from a file
        case "file":

            # get info from comments
            infos = {}
            with open(kwargs["path"], "r") as f:
                for line in f:
                    # parse the comment lines from the file
                    if line.startswith("#"):
                        l = line.split(":")
                        k = l[0]
                        k = k[1:].strip()
                        v = "".join(l[1:])
                        v = v.strip()
                        infos[k] = v
                    else:
                        break
            # retrieve important infos and supplement missing
            if "name" in infos.keys():
                name = infos["name"]
            else:
                warnings.warn(
                    "No name found for externally loaded medium. Setting random one."
                )
                name = "noname_" + "".join(
                    random.choices(string.ascii_uppercase + string.digits, k=6)
                )

            if "description" in infos.keys():
                description = infos["description"]
            else:
                description = None

            if "reference" in infos.keys():
                reference = infos["reference"]
            else:
                reference = None

            # get substances
            substances = read_substances_from_file(kwargs["path"])

            # construct medium
            new_medium = Medium(
                name=name,
                description=description,
                substance_table=substances,
                doi=reference,
            )

            return new_medium

        # unknown case, raise error
        case _:
            raise ValueError(f"Unknown input for parameter how: {how}")


def extract_medium_info_from_model_bigg(row, model: cobra.Model) -> pd.Series:
    """Helper function for :py:func:`read_from_cobra_model`.
    Extracts more information about the medium.

    Args:
        - row (pd.Series):
            A row of the datatable of :py:func:`read_from_cobra_model`.
        - model (cobra.Model):
            The cobra Model

    Returns:
        pd.Series:
            One row for the substance table.
    """

    db_id = row["db_id_EX"].replace("EX_", "").replace("_e", "")
    meta = model.metabolites.get_by_id(db_id + "_e")
    row["name"] = meta.name
    row["formula"] = meta.formula
    row["db_id"] = db_id

    return row


def read_from_cobra_model(
    model: cobra.Model, namespace: Literal["BiGG"] = "BiGG"
) -> Medium:
    """Read and import a medium from a cobra model into a Medium object.

    Args:
        - model (cobra.Model):
            An open cobra Model.

    Returns:
        Medium:
            The imported medium.
    """

    # retrieve the medium from the model
    cobra_medium = model.medium.copy()
    substances = pd.DataFrame(cobra_medium.items(), columns=["db_id_EX", "flux"])
    substances["db_type"] = namespace

    # retrieve additional information
    substances["source"] = model.id
    substances["name"] = None
    substances["formula"] = None

    # transform exchange reacs into the metabolite
    match namespace:
        case "BiGG":
            substances = substances.apply(
                extract_medium_info_from_model_bigg, model=model, axis=1
            )
            # substances['db_id_meta'] = substances.db_id_EX.apply(lambda x: x.removeprefix('EX_'))
            # substances['db_id'] = substances.db_id_meta.apply(lambda x: x.rsplit('_',maxsplit=1)[0])
        case _:
            raise ValueError(f"Unknown input for namespace: {namespace}")

    # reformat table
    substances.drop(columns=["db_id_EX"], inplace=True)
    sorted_cols = ["name", "formula", "flux", "source", "db_id", "db_type"]
    substances = substances.reindex(columns=sorted_cols)

    # construct the medium
    name = f"medium_{model.id}"
    description = f"Medium imported from model {model.id}"
    imported_medium = Medium(name, substances, description)

    return imported_medium


############################################################################
# functions for adding entries to DB
############################################################################

# utility
# -------


def get_last_idx_table(
    tablename: str, connection: sqlite3.Connection, cursor: sqlite3.Cursor
) -> int:
    """Helper function.
    Retrieves the last row id of a specified table of the database.

    Args:
        - tablename (str):
            The name of the table to retrieve the last row id from
        - connection (sqlite3.Connection):
            Connection to the database.
        - cursor (sqlite3.Cursor):
            Cursor for the database.

    Returns:
        int:
            The last row ID of the specified table.
    """

    match tablename:
        case "medium":
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM medium")
        case "medium2substance":
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM medium2substance")
        case "substance":
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM substance")
        case "substance2db":
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM substance2db")
        case _:
            sys.exit(f"Unknown table name for database: {tablename}")

    last_rowid_fetched = last_rowid_res.fetchone()
    last_rowid = last_rowid_fetched[0]

    return last_rowid


# add a new medium
# ----------------


def enter_substance_row(
    row: pd.Series, connection: sqlite3.Connection, cursor: sqlite3.Cursor
) -> int:
    """Helper function for :py:func:`~refinegems.classes.medium.enter_medium_into_db`.
    Enters a new entry in the medium2substance table.

    Args:
        - row (pd.Series):
            A row of the pd.DataFrame of the :py:func:`~refinegems.classes.medium.enter_medium_into_db` function.
        - connection (sqlite3.Connection):
            Connection to the database.
        - cursor (sqlite3.Cursor):
            Cursor for the database.

    Returns:
        int:
            The substance ID in the database of the substance.
    """

    # check if a perfect match exists in database
    exact_match_res = cursor.execute(
        "SELECT * FROM substance WHERE substance.formula = ? AND substance.name = ? ",
        (row["formula"], row["name"]),
    )
    exact_match = exact_match_res.fetchone()

    if exact_match:

        # get substance id
        substance_id = exact_match[0]

    # if not promt for similar or create new substance
    else:

        candidates = cursor.execute(
            "SELECT * FROM substance WHERE substance.formula = ? OR substance.name LIKE ? ",
            (row["formula"], row["name"]),
        )
        candidate_list = candidates.fetchall()

        # if candidates were found, let user choose if one of those should be used

        substance_id = None
        if len(candidate_list) > 0:

            print(
                f'Following similar entries have been found for the substance {row["name"]}.'
            )
            for i in candidate_list:
                print(i)
            res = input(
                "If one matches your entry, please enter the ID or enter skip/s: \n"
            )
            while True:
                # skip
                if res in ["skip", "s"]:
                    break
                # match found
                elif res.isnumeric() and int(res) in [_[0] for _ in candidate_list]:
                    substance_id = int(res)
                    break
                # unknown input, try again
                else:
                    res = input("Unknown input. Please try again:\n")

        # if no matching entry has been found, insert new entry into DB
        if not substance_id:
            cursor.execute(
                "INSERT INTO substance VALUES(?,?,?)",
                (
                    None,
                    row["name"],
                    row["formula"],
                ),
            )
            connection.commit()
            # retrieve the newly added ID
            res = cursor.execute("SELECT last_insert_rowid() FROM substance")
            substance_id = res.fetchone()[0]

    return substance_id


def enter_m2s_row(
    row: pd.Series,
    medium_id: int,
    connection: sqlite3.Connection,
    cursor: sqlite3.Cursor,
):
    """Helper function for :py:func:`refinegems.classes.medium.enter_medium_into_db`.
    Enters a new entry in the medium2substance table.

    Args:
        - row (pd.Series):
            A row of the pd.DataFrame of the :py:func:`~refinegems.classes.medium.enter_medium_into_db` function.
        - medium_id (int):
            The row id of the medium.
        - connection (sqlite3.Connection):
            Connection to the database.
        - cursor (sqlite3.Cursor):
            Cursor for the database.
    """

    # check if entry already exists in database
    exact_match_res = cursor.execute(
        "SELECT 1 FROM medium2substance WHERE medium2substance.medium_id = ? AND medium2substance.substance_id = ? ",
        (medium_id, row["substance_id"]),
    )
    exact_match = exact_match_res.fetchone()

    # else
    if not exact_match:
        cursor.execute(
            "INSERT INTO medium2substance VALUES(?,?,?,?)",
            (
                medium_id,
                row["substance_id"],
                row["flux"],
                row["source"],
            ),
        )
        connection.commit()
    else:
        warnings.warn(
            f'Medium2substance connection {medium_id} - {row["substance_id"]} already exists, skipped second assignment.'
        )


def enter_s2db_row(
    row: pd.Series, db_type: str, connection: sqlite3.Connection, cursor: sqlite3.Cursor
):
    """Helper function for :py:func:`~refinegems.classes.medium.enter_medium_into_db`.
    Enters a new entry in the substance2db table after checking if it has yet to be added.

    Args:
        - row (pd.Series):
            A row of the pd.DataFrame of the :py:func:`~refinegems.classes.medium.enter_medium_into_db` function.
        - db_type (str):
            Type of database identifier to be added.
        - connection (sqlite3.Connection):
            Connection to the database.
        - cursor (sqlite3.Cursor):
            Cursor for the database.
    """

    # only need to enter if ID exists
    if row[db_type]:
        # check if mapping already exists
        entry_check_res = cursor.execute(
            "SELECT 1 FROM substance2db WHERE substance2db.substance_id = ? AND db_id = ?",
            (
                row["substance_id"],
                row[db_type],
            ),
        )
        entry_check = entry_check_res.fetchone()
        if not entry_check:
            cursor.execute(
                "INSERT INTO substance2db VALUES(?,?,?)",
                (
                    row["substance_id"],
                    row[db_type],
                    db_type,
                ),
            )
            connection.commit()


def enter_medium_into_db(medium: Medium, database: str = PATH_TO_DB):
    """Enter a new medium to an already existing database.

    Args:
        - medium (Medium):
            A medium object to be added to the database.
        - database (str, optional):
            Path to the database. Defaults to the in-build databse.
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # collect starting IDs
    sid_medium = get_last_idx_table("medium", connection, cursor)
    sid_medium2substance = get_last_idx_table("medium2substance", connection, cursor)
    sid_substance = get_last_idx_table("substance", connection, cursor)
    sid_substance2db = get_last_idx_table("substance2db", connection, cursor)

    # try to add a new medium to the database
    try:

        # name
        #   check, if name already in DB and prompt, if that problem occurs
        #   enter into database, if no problem remains
        name_check_res = cursor.execute(
            "SELECT 1 FROM medium WHERE medium.name = ?", (medium.name,)
        )
        name_check = name_check_res.fetchone()
        if not name_check:
            # medium name is unique
            # add medium to DB
            cursor.execute(
                "INSERT INTO medium VALUES(?,?,?,?)",
                (
                    None,
                    medium.name,
                    medium.description,
                    medium.doi,
                ),
            )
        else:
            # medium name already in DB
            check = input(
                "The name {medium.name} is already in the database.\nDo you want to set a new name? (yes/no): "
            )
            if check in ["y", "yes"]:
                # set a new name
                while True:
                    new_name = input("Please enter a medium name:\n")
                    name_check_2 = cursor.execute(
                        "SELECT 1 FROM medium WHERE medium.name = ?", (new_name,)
                    )
                    if name_check_2.fetchone():
                        print("This name already exists in the database.")
                    else:
                        medium.name = new_name
                        break
            elif check in ["n", "no"]:
                # end program when no new name is set
                print("No new name chosen. Ending the program.")
                sys.exit()
            else:
                # Abort program at unknown input
                sys.exit("Unknown input. Aborting.")

        connection.commit()
        res = cursor.execute("SELECT last_insert_rowid() FROM medium")
        medium_id = res.fetchone()[0]

        # substances
        medium.substance_table["substance_id"] = medium.substance_table.apply(
            enter_substance_row, axis=1, args=(connection, cursor)
        )

        # medium2substances
        medium.substance_table.apply(
            enter_m2s_row, axis=1, args=(medium_id, connection, cursor)
        )

        # substance2db
        avail_data = [
            _ for _ in medium.substance_table.columns if _ in ALLOWED_DATABASE_LINKS
        ]
        for current_db in avail_data:
            medium.substance_table.apply(
                enter_s2db_row, axis=1, args=(current_db, connection, cursor)
            )

        print(
            f"Medium {medium.name} with ID {medium_id} has been added to the database."
        )

    # in case of problems, interupt and revert changes
    except:

        print(f"During execution, the following error occured: \n{sys.exc_info()}")
        print("Reverting changes to database...")

        cursor.execute("DELETE FROM medium WHERE rowid > ?", (sid_medium,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?", (sid_medium2substance,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?", (sid_substance,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?", (sid_substance2db,))

        print("Done")

    # in any case close connection to database
    finally:
        connection.close()


# add a new subset
# ----------------


def add_subset_to_db(
    name: str,
    desc: str,
    subs_dict: dict,
    database: str = PATH_TO_DB,
    default_perc: float = 1.0,
) -> None:
    """Add a new subset to the database.

    Args:
        - name (str):
            Name (Abbreviation) of the new subset. Needs to be unique for the databse.
        - desc (str):
            Description of the new subset.
        - subs_dict (dict):
            Dictionary of the names and percentages for the substances to be included
            in the new subsets. The names should be part of the substance table.
        - database (str, optional):
            Which database to connect to.
            Defaults to PATH_TO_DB.
        - default_perc (float, optional):
            Default percentage to set if None is given in the dictionary.
            Defaults to 1.0.
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # add new subset to subset table
    res = cursor.execute("SELECT 1 FROM subset WHERE name = ?", (name,))
    if not res.fetchone():
        cursor.execute(
            "INSERT INTO subset VALUES(?,?,?)",
            (
                None,
                name,
                desc,
            ),
        )
        connection.commit()
        res = cursor.execute("SELECT last_insert_rowid() FROM subset")
        subset_id = res.fetchone()[0]

        # add the substance connections to subset2substance
        for s, p in subs_dict.items():

            # set default per if none given
            if not p:
                p = default_perc
            # get substance ID
            res = cursor.execute("SELECT id FROM substance WHERE name = ?", (s,))
            subs_id = res.fetchone()
            if subs_id:
                subs_id = subs_id[0]
                # enter the entry into subset2substance
                cursor.execute(
                    "INSERT INTO subset2substance VALUES(?,?,?)",
                    (
                        subset_id,
                        subs_id,
                        p,
                    ),
                )
                connection.commit()
            else:
                mes = f"Substance {s} not in database. Not added. Please check your input."
                warnings.warn(mes)

    else:
        mes = f"Subset name {name} already in database. Cannot add.\nPlease choose a different name or delete the old one."
        warnings.warn(mes, category=UserWarning)

    # close connection to database
    connection.close()


# further database curation
# -------------------------

def update_db_entry_single(
    table: str,
    column: str,
    new_value: Any,
    conditions: dict,
    database: str = PATH_TO_DB,
):
    """Update a single database entry.

    Args:
        - table (str):
            Name of the table to update.
        - column (str):
            Name of the Attribute to change.
        - new_value (Any):
            New value to be set.
        - conditions (dict):
            Further conditions.
        - database (str, optional):
            Which database to change. Defaults to PATH_TO_DB.
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # condition -> sql
    condition_str = " AND ".join([f"{x} = {y}" for x, y in conditions.items()])

    # update the entry
    cursor.execute(
        "UPDATE ? SET ? = ? WHERE ?",
        (
            table,
            column,
            new_value,
            condition_str,
        ),
    )

    # save and close
    connection.commit()
    connection.close()


def enter_db_single_entry(
    table: str, columns: list[str], values: list[Any], database: str = PATH_TO_DB
):
    """Enter a single entry into a database.

    Args:
        - table (str):
            Which table to enter information to.
        - columns (list[str]):
            Name of the columns to add information to.
        - values (list[Any]):
            List of new values for the columns.
        - database (str, optional):
            Database to add a row to. Defaults to PATH_TO_DB.
    """
    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # add new entry to a specific table
    columns_str = "(" + ", ".join(columns) + ")"
    values_str = "(" + ", ".join(values) + ")"
    cursor.execute(
        "INSERT INTO ? ? VALUES ?",
        (
            table,
            columns_str,
            values_str,
        ),
    )

    # save and close
    connection.commit()
    connection.close()


def generate_update_query(row: pd.Series) -> str:
    """Helper function for :py:func:`~refinegems.classes.medium.update_db_multi`.
    Generates an update SQL query for the provided table

    Args:
        - row (pd.Series):
            Series containing the row of a DataFrame to be used to update a table in a database
            with columns table | column | new_value | conditions

    Returns:
        str:
            SQL query to be used to update a table in a database with the provided data
    """
    colorama_init(autoreset=True)
    table = row["table"]
    conditions_dict = {
        k: v for k, v in [_.split("=") for _ in row["conditions"].split(";")]
    }  # condition (str) : a=x,b=y,....
    new_value = (
        f'\'{row["new_value"]}\'' if type(row["new_value"]) == str else row["new_value"]
    )

    update_query = f'UPDATE {table} SET {row["column"]} = {new_value} WHERE '

    match table:

        case "medium2substance":
            if all(_ in conditions_dict.keys() for _ in ["medium", "substance"]):
                update_query += f"""medium_id = (SELECT medium.id FROM medium WHERE medium.name = \'{conditions_dict.get("medium")}\') \
                AND substance_id = (SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\')\
                """
            else:
                raise ValueError(
                    f"{Fore.MAGENTA}No medium and/or substance keys specified. Chosen table {table} cannot be updated!"
                )

        case "substance2db":
            if "substance" in conditions_dict.keys():
                update_query += f'substance_id = (SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\')'
            else:
                raise ValueError(
                    f"{Fore.MAGENTA}No substance key specified. Chosen table {table} cannot be updated!"
                )

        case _:
            conditions_str = " AND ".join(row["conditions"].split(";"))
            update_query += f"{conditions_str}"

    return update_query


def generate_insert_query(row: pd.Series, cursor) -> str:
    """Helper function for :py:func:`~refinegems.classes.medium.update_db_multi`. Generate the SQL string
    for inserting a new line into the database based on a row of the table.

    Args:
        - row (pd.Series):
            One row of the table of the parent function.

    Returns:
        str:
            The constructed SQL string.
    """
    colorama_init(autoreset=True)
    table = row["table"]

    # condition (str) : a=x,b=y,....
    if row["conditions"] and row["conditions"] not in ["", "-"]:
        conditions_dict = {
            k: v for k, v in [_.split("=") for _ in row["conditions"].split(";")]
        }
    else:
        conditions_dict = None

    # Gather column & new_value inputs in SQL-readable format
    columns_str = f'({row["column"]})'
    value_str = row["new_value"]  # Assign new_value to a variable
    if type(value_str) == str:  # Check if new_value is a string
        if "," in value_str:  # Check if new_value ist actually a list
            # Get all string values as string
            value_str = ", ".join(
                [
                    "'" + _.strip() + "'" if not FLOAT_REGEX.fullmatch(_) else _
                    for _ in row["new_value"].split(",")
                ]
            )
        else:
            value_str = f"'{ value_str.strip()}'"  # new_value is a string

    insert_query = f'INSERT INTO {row["table"]} {columns_str} VALUES '

    match table:

        case "medium2substance":
            if not conditions_dict:
                raise UnboundLocalError(
                    f"{Fore.MAGENTA}No conditions column was found in the provided DataFrame. Chosen table {table} cannot be updated."
                )
            if all(_ in conditions_dict.keys() for _ in ["medium", "substance"]):
                insert_query += f"""(\
                    (SELECT medium.id FROM medium WHERE medium.name = \'{conditions_dict.get("medium")}\'), \
                    (SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\'), \
                    {value_str}\
                    )\
                """
            else:
                raise ValueError(
                    f"{Fore.MAGENTA}No medium and/or substance keys specified. Chosen table {table} cannot be updated!"
                )

        case "subset2substance":
            if not conditions_dict:
                raise UnboundLocalError(
                    f"{Fore.MAGENTA}No conditions column was found in the provided DataFrame. Chosen table {table} cannot be updated."
                )
            if all(_ in conditions_dict.keys() for _ in ["subset", "substance"]):
                insert_query += f"""(\
                    (SELECT subset.id FROM subset WHERE subset.name = \'{conditions_dict.get("subset")}\'), \
                    (SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\'), \
                    {value_str}\
                    )\
                """
            else:
                raise ValueError(
                    f"{Fore.MAGENTA}No medium and/or substance keys specified. Chosen table {table} cannot be updated!"
                )

        case "substance2db":
            if not conditions_dict:
                raise UnboundLocalError(
                    f"{Fore.MAGENTA}No conditions column was found in the provided DataFrame. Chosen table {table} cannot be updated."
                )
            if "substance" in conditions_dict.keys():
                res = cursor.execute(
                    f"""SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\'"""
                )
                substance_id = res.fetchone()[0]
                ins = (
                    cursor.execute(
                        f"""SELECT * FROM substance2db WHERE substance2db.substance_id = {substance_id} AND substance2db.db_id = \'{row["new_value"].split(",")[0].strip()}\'"""
                    )
                    .fetchone()[2]
                    .split("+")
                )
                if len(ins) > 0:
                    cur = row["new_value"].split(",")[1].strip().split("+")
                    missing = [_ for _ in ins if _ not in cur]
                    if len(missing) > 0:
                        update_value = "+".join(sorted(cur + missing))
                        insert_query = f"""UPDATE substance2db SET db_type = \'{update_value}\' WHERE substance2db.substance_id = {substance_id} AND szbbstance2db.db_id = \'{row["new_value"].split(",")[0].strip()}\'"""
                    else:
                        return None
                # insert can be done without any problems
                else:
                    insert_query += f"""(\
                        (SELECT substance.id FROM substance WHERE substance.name = \'{conditions_dict.get("substance")}\'), \
                        {value_str}\
                        )\
                    """

            else:
                raise ValueError(
                    f"{Fore.MAGENTA}No substance key specified. Chosen table {table} cannot be updated!"
                )

        case _:
            insert_query += f"({value_str})"

    return insert_query


def update_db_multi(
    data: pd.DataFrame, update_entries: bool, database: str = PATH_TO_DB
):
    """Updates/Inserts multiple entries in a table from the specified database.
    Given table should have the format:

      row :  table | column | new_value | conditions

    Notes:

    - multiple columns and values are lists with a "," and no whitespaces
    - conditions are listed like: a=x;b=y;...

        - conditions separated by ';'
        - column and value separated by '='
        - no whitespaces

    Args:
        - data (pd.DataFrame):
            DataFrame containing the columns table | column | new_value | conditions
        - update_entries (bool):
            Boolean to determine whether entries should be inserted or updated. False means insert.
        - database (str, optional):
            Path to a database. Defaults to PATH_TO_DB.
    """
    colorama_init(autoreset=True)

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # iterate over the input table
    for idx, row in data.iterrows():

        # row :  table | column | new_value | conditions

        # if conditions given, update
        if row["conditions"] and row["conditions"] not in ["", "-"] and update_entries:
            query = generate_update_query(row)
        # else, insert new values
        else:
            query = generate_insert_query(row, cursor)
            if not query:
                continue

        # update the entry
        try:
            cursor.execute(query)
        except sqlite3.IntegrityError as ie:
            print(f"{Fore.MAGENTA}{ie}")
            print(
                f'Ocurred with: column={row["column"]}, new_value={row["new_value"]}, condition={row["conditions"]}'
            )
            continue

    # save and close
    connection.commit()
    connection.close()


# Function to extract SQL schema with updated SBO/media tables
def updated_db_to_schema(directory: str = "../data/database", inplace: bool = False):
    """Extracts the SQL schema from the database data.db & Transfers it into an SQL file

    Args:
        - directory(str,optional):
            Path to the directory of the updated DB.
            Defaults to '../data/database'.
        - inplace(bool, optional):
            If True, uses the default sql-file name, otherwise extends it with the prefix ``updated``.
    """

    # Not needed to be included in Schema
    NOT_TO_SCHEMA = [
        "BEGIN TRANSACTION;",
        "COMMIT;",
        "bigg_to_sbo",
        "ec_to_sbo",
        "bigg_metabolites",
        "bigg_reactions",
        "modelseed_compounds",
    ]
    counter = 0  # To count rows in newly generated file

    if inplace:
        filename = "media_db.sql"
    else:
        filename = "updated_media_db.sql"

    conn = sqlite3.connect(PATH_TO_DB)
    with open(Path(directory, filename), "w") as file:
        for line in iterdump(conn):
            if not (any(map(lambda x: x in line, NOT_TO_SCHEMA))):
                if "CREATE TABLE" in line and counter != 0:
                    file.write(f"\n\n{line}\n")
                else:
                    file.write(f"{line}\n")
                counter += 1
    conn.close()


############################################################################
# working with models
############################################################################

def medium_to_model(
    model: cobra.Model,
    medium: Medium,
    namespace: str = "BiGG",
    default_flux: float = 10.0,
    replace: bool = False,
    double_o2: bool = True,
    add: bool = True,
) -> Union[None, dict[str, float]]:
    """Add a medium to a COBRApy model.

    Args:
        - model (cobra.Model):
            A model loaded with COBRApy.
        - medium (Medium):
            A refinegems Medium object.
        - namespace (str, optional):
            String to set the namespace to use for the model IDs. Defaults to 'BiGG'.
        - default_flux (float, optional):
            Set a default flux for NaN values or all. Defaults to 10.0.
        - replace (bool, optional):
            Option to replace existing flux values with the default if set to True. Defaults to False.
        - double_o2 (bool, optional):
            Double the oxygen amount in the medium. Defaults to True.
        - add (bool, optional):
            If True, adds the medium to the model, else the exported medium is returned. Defaults to True.

    Returns:
        Union[None, dict[str, float]]:
            Either none or the exported medium.
    """

    # export medium to cobra
    exported_medium = medium.export_to_cobra(
        namespace=namespace,
        default_flux=default_flux,
        replace=replace,
        double_o2=double_o2,
    )

    # remove exchanges that do not exist in model
    model_exchanges = [_.id for _ in model.exchanges]
    exported_medium = {k: v for k, v in exported_medium.items() if k in model_exchanges}

    if add:
        # add to model
        model.medium = exported_medium
        return
    else:
        return exported_medium




