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

from collections.abc import Iterable
from enum import Enum
from sqlite3 import Error
from tqdm import tqdm
from pathlib import Path
from importlib.resources import files
from typing import Union

from ..warnings import BIGG_TABLES, MEDIA_TABLES, MODELSEED_REQUIRED_COLUMNS, database_table_problems

################################################################################
# setup logging
################################################################################

logger = logging.getLogger(__name__)

################################################################################
# variables
################################################################################

PATH_TO_DB_FOLDER = files("refinegems.data.database")  #: :meta hide-value:
PATH_TO_DB = PATH_TO_DB_FOLDER.joinpath("data.db")  #: :meta hide-value:
VERSION_FILE = PATH_TO_DB_FOLDER.joinpath(
    "current_bigg_db_version.txt"
)  #: :meta hide-value:
VERSION_URL = "http://bigg.ucsd.edu/api/v2/database_version"  #: :meta:
MODELSEED_COMPOUNDS_CONTENTS_URL = (
    "https://api.github.com/repos/ModelSEED/ModelSEEDDatabase/contents/Biochemistry?ref=master"
)


MNX_FTP_BASE_PATH = "https://www.metanetx.org/ftp/latest/"
mnx_db_namespace = {
    "reac_prop": (
        f"{MNX_FTP_BASE_PATH}reac_prop.tsv",
        ["id", "mnx_equation", "reference", "ec-code", "is_balanced", "is_transport"],
    ),
    "reac_xref": (
        f"{MNX_FTP_BASE_PATH}reac_xref.tsv",
        ["source", "id", "description"],
    ),
    "chem_prop": (
        f"{MNX_FTP_BASE_PATH}chem_prop.tsv",
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
        f"{MNX_FTP_BASE_PATH}chem_xref.tsv",
        ["source", "id", "description"],
    ),
}

DATABASE_NAMESPACES = (
    "BiGG",
    "ChEBI",
    "KEGG",
    "MetaCyc",
    "MetaNetX",
    "SEED",
    "VMH",
)  #: :meta:

_NAMESPACE_ALIASES = {
    "bigg": "BiGG",
    "chebi": "ChEBI",
    "kegg": "KEGG",
    "metacyc": "MetaCyc",
    "metanetx": "MetaNetX",
    "seed": "SEED",
    "vmh": "VMH",
} #: :meta:

################################################################################
# functions
################################################################################


def canonicalize_namespace(value: object) -> str:
    """Return a validated canonical database namespace.

    Whitespace and accidentally repeated lines are normalized. Legacy combined
    values must be handled with :func:`expand_namespace_memberships`.
    """
    parts = {part.strip() for part in str(value or "").splitlines() if part.strip()}
    if len(parts) != 1:
        raise ValueError(f"Expected one database namespace, received {value!r}")
    raw = next(iter(parts))
    if "+" in raw:
        raise ValueError(f"Combined database namespace is not canonical: {raw!r}")
    try:
        return _NAMESPACE_ALIASES[raw.casefold()]
    except KeyError as error:
        raise ValueError(f"Unknown database namespace: {raw!r}") from error


def expand_namespace_memberships(value: object) -> tuple[str, ...]:
    """Expand legacy combined or repeated input into canonical memberships."""
    raw_parts = []
    for line in str(value or "").splitlines():
        raw_parts.extend(part for part in line.split("+") if part.strip())
    memberships = {canonicalize_namespace(part) for part in raw_parts}
    if not memberships:
        raise ValueError("Database namespace cannot be empty")
    return tuple(name for name in DATABASE_NAMESPACES if name in memberships)


def namespace_matches(value: object, namespace: str) -> bool:
    """Check exact membership while accepting legacy combined input."""
    canonical = canonicalize_namespace(namespace)
    return canonical in expand_namespace_memberships(value)


def expand_namespace_columns(columns: Iterable[object]) -> tuple[str, ...]:
    """Return canonical database columns from an arbitrary column iterable."""
    found = set()
    for column in columns:
        try:
            found.update(expand_namespace_memberships(column))
        except ValueError:
            continue
    return tuple(name for name in DATABASE_NAMESPACES if name in found)


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
       - 6 tables with names 'medium', 'substance', 'substance2db', 'medium2substance', 'subset' & 'subset2substance'
       - 1 table with name 'modelseed_compounds'

    Args:
       - db_cursor (sqlite3.Cursor):
          Cursor from open connection to the database (data.db)

    Returns:
       int:
          Corresponding to one of the ValidationCodes
    """
    print("Verifying database...")

    problems = database_table_problems(db_cursor.connection)
    bigg_tables_contained = not any(name in problems for name in BIGG_TABLES)
    media_tables_contained = not any(name in problems for name in MEDIA_TABLES)
    modelseed_cmpd_tbl_contained = "modelseed_compounds" not in problems

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


def create_media_database(db_cursor: sqlite3.Cursor, reload: bool = False):
    """Create missing or empty media and subset tables from the bundled media_db.sql.

    Existing non-empty tables and their contents are preserved unless ``reload=True``.
    Reloading replaces all six tables and discards local changes. Changes are
    rolled back if loading fails.

    Args:
       - db_cursor (sqlite3.Cursor):
          Cursor from an open connection to the database (data.db).
       - reload (bool, optional):
          Replace all media and subset tables. Defaults to False.
    """
    tables = MEDIA_TABLES
    problems = database_table_problems(db_cursor.connection)
    selected = tables if reload else tuple(name for name in tables if name in problems)
    if not selected:
        return

    # Load the SQL separately so inserts for existing tables are never replayed.
    # Create all missing schemas before copying rows to satisfy foreign keys.
    source = sqlite3.connect(":memory:")
    try:
        with open(files("refinegems.data.database").joinpath("media_db.sql")) as schema:
            source.executescript(schema.read())
        print("Reloading media tables..." if reload else "Adding missing media tables...")
        db_cursor.execute("SAVEPOINT media_tables")
        try:
            if reload:
                for name in reversed(tables):
                    db_cursor.execute(f'DROP TABLE IF EXISTS "{name}"')
            for name in selected:
                if not reload and problems[name] == "empty":
                    continue
                schema_sql = source.execute(
                    "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", (name,)
                ).fetchone()[0]
                db_cursor.execute(schema_sql)
            for name in selected:
                rows = source.execute(f'SELECT * FROM "{name}"')
                placeholders = ", ".join("?" for _ in rows.description)
                db_cursor.executemany(f'INSERT INTO "{name}" VALUES ({placeholders})', rows)
            db_cursor.execute("RELEASE SAVEPOINT media_tables")
        except Error:
            db_cursor.execute("ROLLBACK TO SAVEPOINT media_tables")
            db_cursor.execute("RELEASE SAVEPOINT media_tables")
            raise
    finally:
        source.close()


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
                value = re.sub(r"^meta:", "", value, re.IGNORECASE)
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
        logger.warning(
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
    """Combines ModelSEED's compound TSV partitions into the compounds table.

    All downloads and their required columns are checked before replacing the
    existing table. HTTP failures and invalid compound data raise an exception.

    Args:
       - db_connection (sqlite3.Connection):
          Open connection to the database (data.db)
    """
    print("Adding the ModelSEED compounds table...")
    # ModelSEED now distributes compounds in numbered TSV partitions. Discover
    # them rather than assuming a fixed number, excluding provenance TSVs.
    response = requests.get(MODELSEED_COMPOUNDS_CONTENTS_URL, timeout=60)
    response.raise_for_status()
    compound_files = sorted(
        (
            entry for entry in response.json()
            if entry["type"] == "file"
            and re.fullmatch(r"compound_\d+\.tsv", entry["name"])
        ),
        key=lambda entry: entry["name"],
    )
    if not compound_files:
        raise ValueError("No ModelSEED compound TSV partitions found.")

    frames = []
    for entry in compound_files:
        response = requests.get(entry["download_url"], timeout=60)
        response.raise_for_status()
        frame = pd.read_csv(
            io.StringIO(response.text), sep="\t", dtype={"linked_compound": str}
        )
        if frame.empty or not MODELSEED_REQUIRED_COLUMNS.issubset(frame.columns):
            raise ValueError(f"Invalid ModelSEED compounds data in {entry['name']}.")
        frames.append(frame)

    modelseed_df = pd.concat(frames, ignore_index=True)
    if modelseed_df["id"].isna().any() or modelseed_df["id"].duplicated().any():
        raise ValueError("ModelSEED compound IDs must be non-null and unique.")
    # Leave the existing table untouched if any download or validation fails.
    modelseed_df.to_sql(
        "modelseed_compounds",
        db_connection,
        if_exists="replace",
        index=False,
        dtype={"id": "TEXT PRIMARY KEY", "mass": "NUMERIC", "linked_compound": "TEXT"},
    )


def initialise_database(reload_media: bool = False):
    """Initialises/updates the database (data.db)

    Missing or empty media and subset tables are populated by default. Set
    ``reload_media=True`` to replace all six tables with the bundled SQL
    contents, discarding local changes.

    After initialisation the database contains:

    - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
    - 6 tables with names 'medium', 'substance', 'medium2substance', 'substance2db', 'subset' & 'subset2substance'
    - 1 table with name 'modelseed_compounds'
    """
    print("Initialising database...")
    con = sqlite3.connect(PATH_TO_DB)
    try:
        create_media_database(con.cursor(), reload=reload_media)
        problems = database_table_problems(con)
        get_latest_bigg_databases(
            con, is_missing=any(name in problems for name in BIGG_TABLES)
        )
        if "modelseed_compounds" in problems:
            get_modelseed_compounds_database(con)
        remaining = database_table_problems(con)
        if remaining:
            raise RuntimeError(
                "Database initialisation incomplete: "
                + ", ".join(f"{reason} {name}" for name, reason in remaining.items())
            )
        print(validation_messages[ValidationCodes.COMPLETE])
    finally:
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


def reset_database(database: Union[Path, str] = PATH_TO_DB, for_release: bool = False):
    """Remove tables for certain databases to allow pushing of the database
    to GitHub (reduce size).

    Args:
        - database (Path | str, optional):
            Path to the database. Defaults to PATH_TO_DB, the in-build database.
        - for_release (bool, optional):
            Reload media/subset tables from SQL, remove all downloadable
            namespaces and compact the file. Discards local media changes.
            Defaults to False (remove only MetaNetX tables).
    """
    database = Path(database)
    con = sqlite3.connect(database)
    try:
        with con:
            cursor = con.cursor()
            # Keep the SQL refresh and table removal in the same transaction.
            cursor.execute("BEGIN")
            if for_release:
                create_media_database(cursor, reload=True)
            tables = [
                "mnx_chem_xref", "mnx_reac_xref", "mnx_chem_prop", "mnx_reac_prop"
            ]
            if for_release:
                tables.extend(["bigg_metabolites", "bigg_reactions", "modelseed_compounds"])
            for name in tables:
                cursor.execute(f'DROP TABLE IF EXISTS "{name}"')
        if for_release:
            con.execute("VACUUM")
    finally:
        con.close()
