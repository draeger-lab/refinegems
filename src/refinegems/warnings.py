"""Lightweight runtime warnings for refineGEMs."""

import sqlite3
import warnings
from contextlib import closing
from pathlib import Path


def warn_if_database_uninitialised(database: Path | None = None):
    """Warn about missing, empty or unreadable required database tables."""
    if database is None:
        database = Path(Path(__file__).parent, "data", "database", "data.db")
    database = Path(database)
    tables = (
        "medium", "substance", "substance2db", "medium2substance",
        "subset", "subset2substance", "bigg_metabolites", "bigg_reactions",
        "modelseed_compounds",
    )
    problems = []
    try:
        # mode=ro also prevents SQLite from creating a missing database.
        with closing(sqlite3.connect(f"{database.resolve().as_uri()}?mode=ro", uri=True)) as connection:
            existing = {
                row[0] for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            for name in tables:
                if name not in existing:
                    problems.append(name)
                elif connection.execute(f'SELECT 1 FROM "{name}" LIMIT 1').fetchone() is None:
                    problems.append(name)
            if "modelseed_compounds" in existing:
                columns = {
                    row[1] for row in connection.execute("PRAGMA table_info(modelseed_compounds)")
                }
                if not {"id", "aliases", "formula", "charge"}.issubset(columns):
                    problems.append("invalid modelseed_compounds schema")
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
