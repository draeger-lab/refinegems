#!/usr/bin/env python
""" Provides functions to load a modelpath with libSBML and cobrapy

Depending on the application the model needs to be loaded with cobra (memote)
or with libSBML (activation of groups).
"""

import sys
import cobra
import os
import memote
import json
from libsbml import *# SBMLReader, GroupsExtension, FbcModelPlugin, CVTerm, writeSBML
from Bio import SeqIO
from bioservices import KEGG
from datetime import datetime
from bs4 import BeautifulSoup

__author__ = "Famke Baeuerle"

def load_model_cobra(modelpath):
    """loads model using cobrapy

    Args:
        modelpath (Str): Path to GEM

    Returns:
        cobra-model: loaded model by cobrapy
    """
    mod = cobra.io.read_sbml_model(modelpath)
    return mod


def load_model_libsbml(modelpath):
    """loads model using libsbml

    Args:
        modelpath (Str): Path to GEM

    Returns:
        libsbml-model: loaded model by libsbml
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath) #read from file
    mod = read.getModel()
    return mod

def load_document_libsbml(modelpath):
    """loads model document using libsbml

    Args:
        modelpath (Str): Path to GEM

    Returns:
        libsbml-document: loaded document by libsbml
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath) #read from file
    return read