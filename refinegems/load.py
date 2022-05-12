#!/usr/bin/env python
""" Provides functions to load a modelpath with libSBML and cobrapy

Depending on the application the model needs to be loaded with cobra (memote)
or with libSBML (activation of groups).
"""

import cobra
import os
import pandas as pd
from libsbml import SBMLReader, writeSBMLToFile

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
    read = reader.readSBMLFromFile(modelpath)  # read from file
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
    read = reader.readSBMLFromFile(modelpath)  # read from file
    return read


def load_medium_custom(mediumpath):
    """Helper function to read medium csv

    Args:
        mediumpath (Str): path to csv file with medium

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium['BiGG_R'] = 'R_EX_' + medium['BiGG'] + '_e'
    medium['BiGG_EX'] = 'EX_' + medium['BiGG'] + '_e'
    return medium


def load_medium_from_db(mediumpath, mediumname):
    """Helper function to read standard media_db.csv

    Args:
        mediumpath (Str): path to csv file with medium database
        mediumname (Str): name of medium to test growth on

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium = medium.loc[medium['medium'] == mediumname]
    medium['BiGG_R'] = 'R_EX_' + medium['BiGG'] + '_e'
    medium['BiGG_EX'] = 'EX_' + medium['BiGG'] + '_e'
    return medium


def load_all_media_from_db(mediumpath):
    """Helper function to extract media definitions from media_db.csv

    Args:
        mediumpath (Str): path to csv file with medium database

    Returns:
        df: pandas dataframe of csv with metabs added as BiGG_EX exchange reactions
    """
    media = pd.read_csv(mediumpath, sep=';')
    media['BiGG_R'] = 'R_EX_' + media['BiGG'] + '_e'
    media['BiGG_EX'] = 'EX_' + media['BiGG'] + '_e'

    media['group'] = media['medium'].ne(media['medium'].shift()).cumsum()
    grouped = media.groupby('group')
    media_dfs = []
    for name, data in grouped:
        media_dfs.append(data.reset_index(drop=True))
    return media_dfs


def write_to_file(model, new_filename):
    """Writes modified model to new file

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename / path for modified model
    """
    new_document = model.getSBMLDocument()
    writeSBMLToFile(new_document, new_filename)
    print("Modified model written to " + new_filename)


def write_report(dataframe, filepath):
    """Writes reports stored in dataframes to xlsx file

    Args:
        dataframe (pd.DataFrame): table containing output
        filepath (string): path to file with filename
    """
    writer = pd.ExcelWriter(str(os.path.abspath('.')) + '/' + filepath)
    dataframe.to_excel(writer)
    writer.save()
