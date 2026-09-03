#!/usr/bin/env python
"""This module contains decorator functions."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import logging
import warnings

from functools import wraps, partial

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################


def template(func):
    """A decorator for template functions.

    Template functions are not executable and are mainly for giving developers
    a sense of how functions for the functionality this template is for should look like.
    """

    @wraps(func)
    def wrapper():
        raise RuntimeError(
            "This function is a template for developers.\nIt cannot be executed."
        )

    return wrapper


def implement(func):
    """A decorator for functions that need to be implemented.

    Used to give a hint to other developers, that this function is not yet implemented,
    but should be to make the program executable.
    """

    @wraps(func)
    def wrapper():
        raise NotImplementedError(
            "The current function is just a placeholder and will be implement in the fucture."
        )

    doc_extension = "\n\n.. note::\n\n   *Will be coming in a future release.*"
    if wrapper.__doc__:
        wrapper.__doc__ += doc_extension
    else:
        wrapper.__doc__ = doc_extension
    return wrapper


def debug(func):
    """A decorator for debug functions.

    Debug functions are for giving developers detailed output for debugging the code.
    """

    @wraps(func)
    def wrapper():
        raise RuntimeError("This function is for debugging for the developers.")

    return wrapper


def deprecate(func=None, note: str = None):
    """A decorator to tell the user, that the function will soon be deprecated.

    Used to give hints to users, that their code will not work as expected, if they update
    to a newer version.
    """

    if func is None:
        return partial(deprecate, note=note)

    @wraps(func)
    def wrapper(*args, **kwargs):
        mes = f"This function will be deprecated in the next major update: {func.__name__}"
        if note:
            mes = mes + "\n" + f"Additional notes: {note}"
        logging.warning(mes)
        return func(*args, **kwargs)

    return wrapper


# with help from Copilot 
def suppress_log_message(logger_name, level, message_substring):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            class MessageFilter(logging.Filter):
                def filter(self, record):
                    return not (record.levelno == level and message_substring in record.getMessage())
            filt = MessageFilter()
            logger.addFilter(filt)
            try:
                return func(*args, **kwargs)
            finally:
                logger.removeFilter(filt)
        return wrapper
    return decorator

# written with help from Copilot
def suppress_warning(msg_substring):
    """A decorator to suppress warnings that contain a specific substring.

    Args:
        - msg_substring (str): 
            The substring to look for in warning messages.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings():
                def custom_filter(message, category, filename, lineno, file=None, line=None):
                    if msg_substring in str(message):
                        return
                    return warnings._showwarnmsg_impl(warnings.WarningMessage(
                        message, category, filename, lineno, file, line))
                warnings.showwarning = custom_filter
                return func(*args, **kwargs)
        return wrapper
    return decorator
