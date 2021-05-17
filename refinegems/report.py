#!/usr/bin/env python
""" Provides functions to write properties of a GEM into a report

Initial analysis covers number of reactions, metabolites and genes. 
This part of the script also adds the memote score to the script 
when given to the functions.
"""

import sys
import click
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

def write_report(name, title, react, metab, genes, memote):
    """short report with key features of GEM
        filename 'name_report_day'

    Args:
        name (Str): name of the GEM
        title (Str): model with model Id
        react (Str): number of reactions
        metab (Str): number of metabolites
        genes (Str): number of genes
        memote (Str): memote score
    """
    now = datetime.now()
    day = now.strftime("%Y%m%d")
    filename = str(name) + '_report_' + str(day)
    original_stdout = sys.stdout  # Save a reference to the original standard output
    with open(filename, 'w') as f:
        sys.stdout = f  # Change the standard output to the file we created.
        print(title)
        print('--- report created on ' + 
                str(now.strftime("%Y/%m/%d, %H:%M")) + '---')
        print(react)
        print(metab)
        print(genes)
        print(memote)
        sys.stdout = original_stdout  # Reset the standard output to its original value
    print('Wrote report to ' + filename)

        name = model.getId()
    title = 'This is a model of ' + model.getId()
    reactions = '%i reactions' % model.getNumReactions()
    metabolites = '%i metabolites' % model.getNumSpecies()
    genes = '%i genes' % len(model.getPlugin(0).getListOfGeneProducts())


def write_to_cl(name, title, react, metab, genes, memote):
    """wrire key features to command line

    Args:
        name (Str): name of the GEM
        title (Str): model with model Id
        react (Str): number of reactions
        metab (Str): number of metabolites
        genes (Str): number of genes
        memote (Str): memote score
    """
    now = datetime.now()
    print(title)
    print('--- report created on ' +
            str(now.strftime("%Y/%m/%d, %H:%M")) + ' ---')
    print(react)
    print(metab)
    print(genes)
    print(memote)