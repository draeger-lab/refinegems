"""Logging configuration for interactive refineGEMs use and the CLI."""

################################################################################
# requirements
################################################################################

import logging
from copy import copy

import click

################################################################################
# variables
################################################################################

_LEVEL_COLOURS = {
    logging.DEBUG: "bright_blue",
    logging.INFO: "green",
    logging.WARNING: "yellow",
    logging.ERROR: "red",
    logging.CRITICAL: "red",
}
_CONSOLE_HANDLER_NAME = "refinegems.console"

################################################################################
# classes
################################################################################


class LevelFormatter(logging.Formatter):
    """Format diagnostics with optional colour confined to the level keyword."""

    def __init__(self, colour: bool = True):
        """Initialise the diagnostic formatter.

        Args:
            - colour (bool, optional):
                Colour the level keyword. Defaults to True.
        """
        super().__init__("%(levelname)s [%(name)s]: %(message)s")
        self.colour = colour

    def format(self, record: logging.LogRecord) -> str:
        """Format a log record without changing it for other handlers.

        Args:
            - record (logging.LogRecord):
                Diagnostic record to format.

        Returns:
            str:
                Formatted diagnostic, including any exception traceback.
        """
        # Other handlers must receive the original, uncoloured level name.
        record = copy(record)
        if self.colour and record.levelno in _LEVEL_COLOURS:
            record.levelname = click.style(
                record.levelname,
                fg=_LEVEL_COLOURS[record.levelno],
                bold=True if record.levelno == logging.CRITICAL else None,
            )
        return super().format(record)


################################################################################
# functions
################################################################################


def configure_logging(level: str = "INFO", colour: bool = True) -> logging.Logger:
    """Configure refineGEMs console logging, including in Jupyter notebooks.

    Args:
        - level (str, optional):
            Minimum diagnostic level. Defaults to 'INFO'.
        - colour (bool, optional):
            Colour only the level keyword. Defaults to True.
            Set to False for plain output.

    Raises:
        - ValueError: Unknown logging level

    Returns:
        logging.Logger:
            Configured parent logger for refineGEMs modules.
    """
    # Validate before changing an existing logging configuration.
    handler = logging.StreamHandler()
    handler.setLevel(level.upper())
    handler.setFormatter(LevelFormatter(colour=colour))
    handler.set_name(_CONSOLE_HANDLER_NAME)
    logger = logging.getLogger("refinegems")

    # Re-running notebook cells replaces our handler rather than duplicating it.
    for existing in list(logger.handlers):
        if existing.get_name() == _CONSOLE_HANDLER_NAME:
            logger.removeHandler(existing)
            existing.close()
    logger.addHandler(handler)
    logger.setLevel(level.upper())
    logger.propagate = False  # Avoid repeating diagnostics through root handlers.
    return logger
