"""Lightweight runtime warnings for refineGEMs."""

import sqlite3
import warnings
from contextlib import closing
from pathlib import Path


MEDIA_TABLES = (
    "medium", "substance", "subset", "medium2substance",
    "substance2db", "subset2substance",
)
BIGG_TABLES = ("bigg_metabolites", "bigg_reactions")
MODELSEED_REQUIRED_COLUMNS = {
    "id", "name", "formula", "charge", "mass", "linked_compound", "aliases"
}


def database_table_problems(connection: sqlite3.Connection) -> dict[str, str]:
    """Return missing, empty or invalid required tables for warnings and repair."""
    existing = {
        row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    problems = {}
    for name in (*MEDIA_TABLES, *BIGG_TABLES, "modelseed_compounds"):
        if name not in existing:
            problems[name] = "missing"
        elif connection.execute(f'SELECT 1 FROM "{name}" LIMIT 1').fetchone() is None:
            problems[name] = "empty"
    if "modelseed_compounds" in existing:
        columns = {
            row[1] for row in connection.execute("PRAGMA table_info(modelseed_compounds)")
        }
        if not MODELSEED_REQUIRED_COLUMNS.issubset(columns):
            problems["modelseed_compounds"] = "invalid schema"
    return problems


def warn_if_database_uninitialised(database: Path | None = None):
    """Warn about missing, empty or unreadable required database tables."""
    if database is None:
        database = Path(Path(__file__).parent, "data", "database", "data.db")
    database = Path(database)
    try:
        # mode=ro also prevents SQLite from creating a missing database.
        with closing(sqlite3.connect(f"{database.resolve().as_uri()}?mode=ro", uri=True)) as connection:
            problems = [
                f"{reason} {name}" for name, reason in database_table_problems(connection).items()
            ]
    except (sqlite3.Error, OSError, ValueError):
        problems = ["database missing or unreadable"]
    if problems:
        warnings.warn(
            "The refineGEMs database is not fully initialised ("
            + ", ".join(problems)
            + "). Run `refinegems database initialise` before using database-dependent features.",
            RuntimeWarning,
            stacklevel=2,
        )
