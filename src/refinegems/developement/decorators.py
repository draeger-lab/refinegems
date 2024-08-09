#!/usr/bin/env python
"""This module contains decorator functions.
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

def template(func):
    """A decorator for template functions.
    
    Template functions are not executable and are only for developers.
    """
    def wrapper():
        raise RuntimeError('This function is a template for developers.\nIt cannot be executed.')
    return wrapper


def implement(func):
    """A decorator for functions, that should be implemented.
    
    Used to give a hint to other developers, that this function is not yet implemented, 
    but should be to make the program executable.
    """
    def wrapper():
        raise NotImplementedError('The current function is just a placeholder and will be implement in the fucture.')
    return wrapper