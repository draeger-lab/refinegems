#!/usr/bin/env python
"""This module contains decorator functions.
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import warnings
from functools import wraps

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# @TODO: Add note message for docs
def template(func):
    """A decorator for template functions.
    
    Template functions are not executable and are mainly for giving developers 
    a sense of how functions for the functionality this template is for should look like.
    """
    @wraps(func)
    def wrapper():
        raise RuntimeError('This function is a template for developers.\nIt cannot be executed.')
    return wrapper


def implement(func):
    """A decorator for functions that need to be implemented.
    
    Used to give a hint to other developers, that this function is not yet implemented, 
    but should be to make the program executable.
    """
    @wraps(func)
    def wrapper():
        raise NotImplementedError('The current function is just a placeholder and will be implement in the fucture.')
    
    doc_extension = """\n.. note:: \n\t*Will be coming in a future release.*"""
    if wrapper.__doc__:
        wrapper.__doc__ += doc_extension
    else:
        wrapper.__doc__ = doc_extension
    return wrapper
    

# @TODO add an option for alternative / new function
def deprecate(func):
    """A decorator to tell the user, that the function will soon be deprecated.
    
    Used to give hints to users, that their code will not work as expected, if they update 
    to a newer version.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        mes = f'This function will be deprecated in the next major update: {func.__name__}'
        warnings.warn(mes, DeprecationWarning)
        func(*args, **kwargs)
    return wrapper