#!/usr/bin/env python
"""This module contains decorator functions.
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import warnings

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

def template(func):
    """A decorator for template functions.
    
    Template functions are not executable and are mainly for giving developers 
    a sense of how function for the functionality this template is for should look like.
    """
    def wrapper():
        raise RuntimeError('This function is a template for developers.\nIt cannot be executed.')
    return wrapper


def implement(func):
    """A decorator for functions that need to be implemented.
    
    Used to give a hint to other developers, that this function is not yet implemented, 
    but should be to make the program executable.
    """
    def wrapper():
        raise NotImplementedError('The current function is just a placeholder and will be implement in the fucture.')
    return wrapper


# @TODO add an option for alternative / new function
def deprecate(func):
    """A decorator to tell the user, that the function will soon be deprecated.
    
    Used to give hints to users, that their code will not work as expected, if they update 
    to a newer version.
    """
    
    def wrapper(*args, **kwargs):
        mes = f'This function will be deprecated in the next major update: {func.__name__}'
        warnings.warn(mes, DeprecationWarning)
        func(*args, **kwargs)
    return wrapper