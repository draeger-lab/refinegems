#!/usr/bin/env python
"""Variables, functions and more for the developement, extension and maintainance of
the in-build database.

.. note::

   Some functionalities for handling and dealing with the media are in the
   :py:mod:`~refinegems.classes.medium` module.

.. hint::

   Further functions for accessing the database can be found in the :py:mod:`~refinegems.utility.io`
   module, e.g. :py:func:`~refinegems.utility.io.load_a_table_from_database`
"""

__author__ = "Gwendolyn O. Döbel und Carolin Brune"

################################################################################
# requirements
################################################################################

import io
import re
import sqlite3
import requests
import logging
import pandas as pd

from enum import Enum
from sqlite3 import Error
from tqdm import tqdm
from pathlib import Path
from importlib.resources import files
from typing import Union

################################################################################
# variables
################################################################################

PATH_TO_DB_FOLDER = files("refinegems.data.database")  #: :meta hide-value:
PATH_TO_DB = PATH_TO_DB_FOLDER.joinpath("data.db")  #: :meta hide-value:
VERSION_FILE = PATH_TO_DB_FOLDER.joinpath(
    "current_bigg_db_version.txt"
)  #: :meta hide-value:
VERSION_URL = "http://bigg.ucsd.edu/api/v2/database_version"  #: :meta:

mnx_db_namespace = {
    "reac_prop": (
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",
        ["id", "mnx_equation", "reference", "ec-code", "is_balanced", "is_transport"],
    ),
    "reac_xref": (
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        ["source", "id", "description"],
    ),
    "chem_prop": (
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",
        [
            "id",
            "name",
            "reference",
            "formula",
            "charge",
            "mass",
            "InChI",
            "InChIKey",
            "SMILES",
        ],
    ),
    "chem_xref": (
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
        ["source", "id", "description"],
    ),
}

################################################################################
# functions
################################################################################


class ValidationCodes(Enum):
    """Validation codes for the database

    Args:
       - Enum (Enum):
          Provided as input to get a number mapping for the codes
    """

    COMPLETE = (0,)  # All tables are in data.db
    EMPTY = (1,)  # data.db is either empty or incorrect
    BIGG = (2,)  # Only BiGG tables are in data.db
    MEDIA = (3,)  # Only Media tables are in data.db 
    BIGG_MEDIA = (4,)  # Only BiGG and media tables are in data.db 
    MODELSEED_COMPOUNDS = (5,)  # Only ModelSEED compounds table is in data.db
    BIGG_MSEED_COMPOUNDS = (
        6,
    )  # Only Bigg and ModelSEED compounds tables are in data.db
    MEDIA_MSEED_COMPOUNDS = (
        7  # Only media and ModelSEED compounds tables are in data.db # 
    )


validation_messages = {
    ValidationCodes.COMPLETE: "All tables in data up-to-date. Initialisation complete.",
    ValidationCodes.EMPTY: "No table in data. An error must have occurred during initialisation.",
    ValidationCodes.BIGG: "Data only contains the BiGG tables. Please check the remaining tables.",
    ValidationCodes.MEDIA: "Data only contains the media tables. Please check the BiGG and ModelSEED tables.", 
    ValidationCodes.BIGG_MEDIA: "Data only contains the BiGG and media tables. Please check the ModelSEED table.", 
    ValidationCodes.MODELSEED_COMPOUNDS: "Data only contains the ModelSEED table. Please check the BiGG and media tables.", 
    ValidationCodes.BIGG_MSEED_COMPOUNDS: "Data only contains the BiGG and ModelSEED tables. Please check the media tables.",
    ValidationCodes.MEDIA_MSEED_COMPOUNDS: "Data only contains the media and ModelSEED tables. Please check the BiGG tables.", 
}


def is_valid_database(db_cursor: sqlite3.Cursor) -> int:
    """Verifies if database has:

       - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
       - 6 tables with names 'medium', 'substance', 'substance2db' & 'medium2substance', 'subset' & 'subset2substance'
       - 1 table with name 'modelseed_compounds'

    Args:
       - db_cursor (sqlite3.Cursor):
          Cursor from open connection to the database (data.db)

    Returns:
       int:
          Corresponding to one of the ValidationCodes
    """
    print("Verifying database...")

    # Fetches the table names as string tuples from the connected database
    db_cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [string[0] for string in db_cursor.fetchall()]

    bigg_tables_contained = (
        len([s for s in tables if re.match(r"^bigg_(?!to)(.*?)", s, re.IGNORECASE)])
        == 2
    )
    media_tables_contained = (
        len(
            [
                s
                for s in tables
                if re.match(
                    r"^medium(.*?)|^substance(.*?)|^subset(.*?)", s, re.IGNORECASE
                )
            ]
        )
        == 6
    )
    modelseed_cmpd_tbl_contained = (
        len([s for s in tables if s == "modelseed_compounds"]) == 1
    )

    bigg_media_tbls_contained = bigg_tables_contained and media_tables_contained 
    bigg_modelseed_cmpd_tbls_contained = (
        bigg_tables_contained and modelseed_cmpd_tbl_contained
    )
    media_modelseed_cmpd_tbls_contained = (
        media_tables_contained and modelseed_cmpd_tbl_contained 
    )
    all_tables_contained = (
        bigg_media_tbls_contained and modelseed_cmpd_tbl_contained
    )

    if all_tables_contained:
        return ValidationCodes.COMPLETE
    elif bigg_modelseed_cmpd_tbls_contained:
        return ValidationCodes.BIGG_MSEED_COMPOUNDS
    elif media_modelseed_cmpd_tbls_contained:
        return ValidationCodes.MEDIA_MSEED_COMPOUNDS
    elif bigg_media_tbls_contained:
        return ValidationCodes.BIGG_MEDIA 
    elif bigg_tables_contained:
        return ValidationCodes.BIGG
    elif media_tables_contained:
        return ValidationCodes.MEDIA 
    elif modelseed_cmpd_tbl_contained:
        return ValidationCodes.MODELSEED_COMPOUNDS
    else:
        return ValidationCodes.EMPTY


def create_media_database(db_cursor: sqlite3.Cursor):
    """Creates the media database with 4 tables 
    ('medium', 'substance', 'substance2db', 'medium2substance') 
    from file './data/database/media_db.sql'

    Args:
       - db_cursor (sqlite3.Cursor):
          Cursor from open connection to the database (data.db)
    """
    db_cursor.executescript(
        """
                           DROP TABLE IF EXISTS susbtance;
                           DROP TABLE IF EXISTS substance2db;
                           DROP TABLE IF EXISTS medium;
                           DROP TABLE IF EXISTS medium2substance;
                           """
    )

    print("Adding media tables...")
    with open(files("refinegems.data.database").joinpath("media_db.sql")) as schema:
        db_cursor.executescript(schema.read())


def update_bigg_db(latest_version: str, db_connection: sqlite3.Connection) -> dict:
    """Updates the BiGG tables 'bigg_metabolites' & 'bigg_reactions' within a database (data.db)

    Args:
      - latest_version (str):
         String containing the Path to a file with the latest version of the BiGG database
      - db_connection (sqlite3.Connection):
         Open connection to the database (data.db)
    """

    def get_database_links_info_per_row(row: pd.Series):
        """For a single row of the dataframe, extract the database identifier from the collection of database links.

        Args:
            - ow (pd.Series):
               The database links for a single row of the dataframe.

        Returns:
            dict:
               A dictionary with the databses as keys and the IDs as values.
        """

        database_ids = {}
        if isinstance(row, str):
            links = row.split(";")
            for link in links:
                key, value = link.split(":", 1)
                key = key.strip()
                value = value.rsplit("/", 1)[1].strip()
                value = re.sub(r"^(?i)meta:", "", value)
                if key in database_ids.keys():
                    database_ids[key].append(value)
                else:
                    database_ids[key] = [value]

        for k in database_ids.keys():
            database_ids[k] = ", ".join(database_ids[k])

        return database_ids

    def get_database_links(data: pd.DataFrame) -> pd.DataFrame:
        """For the dataframe, extract the database IDs from the database_links column
        into separate column containing the IDs and not the links.

        Args:
            - data (pd.DataFrame):
               The input dataframe.

        Returns:
            pd.DataFrame:
               The edited dataframe
        """

        data = data.join(
            pd.DataFrame(
                [get_database_links_info_per_row(row) for row in data["database_links"]]
            )
        )
        data.drop(columns=["database_links", "model_list"], axis=1, inplace=True)

        return data

    print("Adding BiGG tables...")

    # Store currently used version
    with open(VERSION_FILE, "w") as file:
        file.write(latest_version)

    # Create BiGG metabolites table
    BIGG_MODELS_METABS_URL = (
        "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
    )
    bigg_models_metabs = requests.get(BIGG_MODELS_METABS_URL).text
    bigg_models_metabs_df = pd.read_csv(
        io.StringIO(bigg_models_metabs), dtype=str, sep="\t"
    )
    bigg_models_metabs_df.rename(columns={"bigg_id": "id"}, inplace=True)

    bigg_models_metabs_df = get_database_links(bigg_models_metabs_df)

    bigg_id_duplicates = bigg_models_metabs_df.duplicated(subset=["id"], keep=False)
    bigg_id_duplicates_df = bigg_models_metabs_df[bigg_id_duplicates]
    bigg_id_duplicates_set = set(bigg_id_duplicates_df["id"].tolist())

    bigg_models_metabs_df[~bigg_id_duplicates].to_sql(
        "bigg_metabolites",
        db_connection,
        if_exists="replace",
        index=False,
        dtype={"id": "TEXT PRIMARY KEY"},
    )

    if bigg_id_duplicates_set:
        logging.warning(
            "The BiGG metabolite table contains the following "
            f"{len(bigg_id_duplicates_set)} duplicate(s):\n"
            f"{bigg_id_duplicates_set}\n"
            "Duplicate(s) are completely removed from the table."
        )

    # Create BiGG reactions table
    BIGG_MODELS_REACS_URL = (
        "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt"
    )
    bigg_models_reacs = requests.get(BIGG_MODELS_REACS_URL).text
    bigg_models_reacs_df = pd.read_csv(
        io.StringIO(bigg_models_reacs), dtype=str, sep="\t"
    )
    bigg_models_reacs_df.rename(columns={"bigg_id": "id"}, inplace=True)

    bigg_models_reacs_df = get_database_links(bigg_models_reacs_df)

    bigg_models_reacs_df.to_sql(
        "bigg_reactions",
        db_connection,
        if_exists="replace",
        index=False,
        dtype={"id": "TEXT PRIMARY KEY"},
    )


def get_latest_bigg_databases(
    db_connection: sqlite3.Connection, is_missing: bool = True
):
    """Gets the latest BiGG tables for metabolites & reactions if:

       - No version file is locally available
       - The version in the local version file is NOT the latest
       - No BiGG tables currently exist in the database

    Args:
       - db_connection (sqlite3.Connection):
          Open connection to the database (data.db)
       - is_missing (bool, optional):
          True if no BiGG tables are in the database. Defaults to True.
    """
    # Check if BiGG database had an update
    LATEST_VERSION = requests.get(VERSION_URL).json()["bigg_models_version"]

    if not Path.exists(VERSION_FILE) or is_missing:
        update_bigg_db(LATEST_VERSION, db_connection)

    else:
        with open(VERSION_FILE, "r") as file:
            version = file.readline().strip()

        if version != LATEST_VERSION:
            update_bigg_db(LATEST_VERSION, db_connection)


def get_modelseed_compounds_database(db_connection: sqlite3.Connection):
    """Retrieves the compounds table from ModelSEED from the respective GitHub repository

    Args:
       - db_connection (sqlite3.Connection):
          Open connection to the database (data.db)
    """
    print("Adding the ModelSEED compounds table...")
    MODELSEED_COMPOUNDS_URL = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv"
    modelseed_compounds = requests.get(MODELSEED_COMPOUNDS_URL).text
    modelseed_df = pd.read_csv(
        io.StringIO(modelseed_compounds), sep="\t", dtype={"linked_compound": str}
    )  #'mass': float,
    modelseed_df.to_sql(
        "modelseed_compounds",
        db_connection,
        if_exists="replace",
        index=False,
        dtype={"id": "TEXT PRIMARY KEY", "mass": "NUMERIC", "linked_compound": "TEXT"},
    )


def initialise_database():
    """Initialises/updates the database (data.db)

    After initialisation the database contains:

    - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
    - 6 tables with names 'medium', 'substance', 'medium2substance', 'substance2db', 'subset' & 'subset2substance'
    - 1 table with name 'modelseed_compounds'
    """
    # Initialise empty connection
    con = None

    print("Initialising database...")

    # Try to open connection & get cursor
    try:
        con = sqlite3.connect(PATH_TO_DB)
        cursor = con.cursor()

        validity_code = is_valid_database(cursor)

        if validity_code == ValidationCodes.BIGG:
            print("Only BiGG tables contained in database.")
            create_media_database(cursor)
            get_modelseed_compounds_database(con)

        elif validity_code == ValidationCodes.MEDIA: 
            print("Only media tables contained in database.")
            get_latest_bigg_databases(con)
            get_modelseed_compounds_database(con)

        elif validity_code == ValidationCodes.MODELSEED_COMPOUNDS:
            print("Only ModelSEED compounds table contained in database.")
            create_media_database(cursor)
            get_latest_bigg_databases(con)

        elif validity_code == ValidationCodes.BIGG_MEDIA: 
            print("Only BiGG and media tables contained in database.")
            get_modelseed_compounds_database(con)

        elif validity_code == ValidationCodes.BIGG_MSEED_COMPOUNDS:
            print("Only BiGG and ModelSEED compounds tables contained in database.")
            create_media_database(cursor)

        elif validity_code == ValidationCodes.MEDIA_MSEED_COMPOUNDS: 
            print(
                "Only media and ModelSEED compounds tables contained in database."
            )
            get_latest_bigg_databases(con)

        elif validity_code == ValidationCodes.EMPTY:
            print(
                "Incorrect or empty database. Initialise database with required tables..."
            )
            create_media_database(cursor)
            get_latest_bigg_databases(con)
            get_modelseed_compounds_database(con)

        elif validity_code == ValidationCodes.COMPLETE:
            print("Verifying if BiGG tables are up-to-date...")
            get_latest_bigg_databases(con, False)

    except Error as e:
        print(e)
    finally:
        if con:
            con.close()

        # Validate initialised database
        con = sqlite3.connect(PATH_TO_DB)
        cursor = con.cursor()
        validity_code = is_valid_database(cursor)
        print(validation_messages.get(validity_code))
        con.close()


# about MetaNetX Information
# --------------------------
# DISCLAIMER:
# Database information from MetaNetX
# distributed under https://creativecommons.org/licenses/by/4.0/
# Citation: MetaNetX/MNXref: unified namespace for metabolites and biochemical reactions in the context of metabolic models
#           Sébastien Moretti, Van Du T Tran, Florence Mehl, Mark Ibberson, Marco Pagni
#           Nucleic Acids Research (2021), 49(D1):D570-D574
def update_mnx_namespaces(db: Union[Path, str] = PATH_TO_DB, chunksize: int = 1):
    """Add or update the MetaNetX namespace to/in a database.

    Args:
       - db (Union[Path,str],optional):
          Path to a database to add the namespace to.
          Defaults to the in-build database.
       - chunksize (int, optional):
          Size of the chunk (in kB) to download at once.
          Defaults to 1.
    """

    con = sqlite3.connect(db)
    for name, values in mnx_db_namespace.items():
        link, colnames = values
        mnx_table = []
        for chunk in tqdm(
            pd.read_csv(
                link, sep="\t", comment="#", names=colnames, chunksize=chunksize * 1024
            ),
            desc=f"Downloading {name}",
            unit="B",
        ):
            mnx_table.append(chunk)

        match name:
            # Reaction property table
            case "reac_prop":
                total_len = sum([len(_) for _ in mnx_table])
                with tqdm(total=total_len, unit="entries", desc="Add to DB") as pbar:
                    for i, chunk in enumerate(mnx_table):
                        if i == 0:
                            exists = "replace"
                        else:
                            exists = "append"
                        chunk.to_sql(
                            "mnx_" + name,
                            con,
                            if_exists=exists,
                            index=False,
                            dtype={"id": "TEXT PRIMARY KEY"},
                        )
                        pbar.update(len(chunk))

            # Reaction cross-reference table
            case "reac_xref":
                cursor = con.cursor()
                cursor.execute("DROP TABLE IF EXISTS mnx_reac_xref")
                empty_table = """ CREATE TABLE mnx_reac_xref (
                                 source TEXT,
                                 id TEXT,
                                 description TEXT,
                                 CONSTRAINT PK_mnx_reac_xref PRIMARY KEY (source,id)
                                 FOREIGN KEY(id) REFERENCES mnx_reac_prop(id)
                           );
                           """
                cursor.execute(empty_table)
                total_len = sum([len(_) for _ in mnx_table])
                with tqdm(total=total_len, unit="entries", desc="Add to DB") as pbar:
                    for i, chunk in enumerate(mnx_table):
                        chunk.to_sql(
                            "mnx_" + name, con, if_exists="append", index=False
                        )
                        pbar.update(len(chunk))

            # Metabolite properties table
            case "chem_prop":
                total_len = sum([len(_) for _ in mnx_table])
                with tqdm(total=total_len, unit="entries", desc="Add to DB") as pbar:
                    for i, chunk in enumerate(mnx_table):
                        if i == 0:
                            exists = "replace"
                        else:
                            exists = "append"
                        chunk.to_sql(
                            "mnx_" + name,
                            con,
                            if_exists=exists,
                            index=False,
                            dtype={"id": "TEXT PRIMARY KEY"},
                        )
                        pbar.update(len(chunk))

            # Metabolite reference table
            case "chem_xref":
                total_len = sum([len(_) for _ in mnx_table])
                cursor = con.cursor()
                cursor.execute("DROP TABLE IF EXISTS mnx_chem_xref")
                empty_table = """ CREATE TABLE mnx_chem_xref (
                                 source TEXT,
                                 id TEXT,
                                 description TEXT,
                                 CONSTRAINT PK_mnx_chem_xref PRIMARY KEY (source,id)
                                 FOREIGN KEY(id) REFERENCES mnx_chem_prop(id)
                           );
                           """
                cursor.execute(empty_table)
                with tqdm(total=total_len, unit="entries", desc="Add to DB") as pbar:
                    for i, chunk in enumerate(mnx_table):
                        chunk.to_sql(
                            "mnx_" + name, con, if_exists="append", index=False
                        )
                        pbar.update(len(chunk))
    con.close()


def reset_database(database: Union[Path, str] = PATH_TO_DB):
    """Remove tables for certain databases to allow pushing of the database
    to GitHub (reduce size).

    Args:
        - database (Path | str, optional):
            Path to the database. Defaults to PATH_TO_DB, the in-build database.
    """
    # establish connection
    con = sqlite3.connect(database)
    cursor = con.cursor()
    # remove MetaNetX tables
    cursor.execute("DROP TABLE IF EXISTS mnx_chem_xref")
    cursor.execute("DROP TABLE IF EXISTS mnx_reac_xref")
    cursor.execute("DROP TABLE IF EXISTS mnx_chem_prop")
    cursor.execute("DROP TABLE IF EXISTS mnx_reac_prop")
    # close connection
    con.close()
