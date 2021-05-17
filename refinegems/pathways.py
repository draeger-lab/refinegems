#!/usr/bin/env python
""" Provides functions for adding KEGG reactions as Group Pathways

If your organism occurs in the KEGG database, extract the KEGG reaction ID from the
annotations of your reactions and identify, in which KEGG pathways this reaction occurs. Add
all KEGG pathways for a reaction then as annotations with the biological qualifier ‘OCCURS_IN’
to the respective reaction.
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
from bs4 import BeautifulSoup

__author__ = "Famke Baeuerle"

def load_model_enable_groups(modelpath):
    """loads model as document using libsbml
        enables groups extension

    Args:
        modelpath (Str): Path to GEM

    Returns:
        libsbml-document: loaded document by libsbml
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath) #read from file
    groupextension = GroupsExtension()
    groupURI = groupextension.getURI(3,1,1)
    read.enablePackage(groupURI, "groups", True) #enable groups extension
    read.setPkgRequired('groups', False) #make groups not required
    model = read.getModel()
    return model

def write_to_file(model, new_filename):
    """Writes modified model to new file

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    new_document = model.getSBMLDocument()
    writeSBMLToFile(new_document, new_filename)
    print("Written to " + new_filename)

def extract_kegg_reactions(model):
    """extract KEGG Ids from reactions

    Args:
        model (libsbml-model): model loaded with libsbml

    Returns:
        dict: reaction Id as key and Kegg Id as value
    """
    list_reac = model.getListOfReactions()
    kegg_reactions = {}

    for reaction in list_reac:
        notes_string = reaction.getNotesString()
        soup = BeautifulSoup(notes_string, 'lxml')
        entries = soup.find_all('p')

        for i in range(len(entries)):
            if 'KEGG' in entries[i].text:
                kegg_reactions[reaction.getId()] = entries[i].text[15:]

    return kegg_reactions

def extract_kegg_pathways(kegg_reactions):
    """finds pathway for KEGG reactions

    Args:
        kegg_reactions (dict): reaction Id as key and Kegg Id as value

    Returns:
        dict: reaction Id as key and Kegg Pathway Id as value
    """
    k = KEGG()
    kegg_pathways = {}

    for reaction in kegg_reactions.keys():
        kegg_reaction = k.get(kegg_reactions[reaction])
        #print(kegg_reaction)
        dbentry = k.parse(kegg_reaction)
        # sometimes parse does not work -> try and except
        try:
            pathways = [x for x in dbentry['PATHWAY']]
        except:
            pathways = []
        kegg_pathways[reaction] = pathways

    return kegg_pathways

def add_kegg_pathways(model, kegg_pathways):
    """add KEGG reactions as BQB_OCCURS_IN

    Args:
        model (libsbml-model): model loaded with libsbml
        kegg_pathways (dict): reaction Id as key and Kegg Pathway Id as value

    Returns:
        libsbml-model: modified model with Kegg pathways
    """
    list_reac = model.getListOfReactions()

    for reaction in list_reac:
        if reaction.getId() in kegg_pathways.keys():
            for path in kegg_pathways[reaction.getId()]:
                url = "https://identifiers.org/kegg.pathway:" + path
                cv = CVTerm()
                cv.setQualifierType(BIOLOGICAL_QUALIFIER)
                cv.setBiologicalQualifierType(BQB_OCCURS_IN)
                cv.addResource(url)
                reaction.addCVTerm(cv)
    
    # works but better write this somewhere else?
    return model

def get_pathway_groups(kegg_pathways):
    """group reaction into pathways

    Args:
        kegg_pathways (dict): reaction Id as key and Kegg Pathway Id as value

    Returns:
        dict: Kegg Pathway Id as key and reactions Ids as values
    """
    pathway_groups = {}
    for reaction in kegg_pathways.keys():
        for path in kegg_pathways[reaction]:
            if path not in pathway_groups.keys():
                pathway_groups[path] = [reaction]
            else:
                pathway_groups[path].append(reaction)
    return pathway_groups

def create_pathway_groups(model, pathway_groups):
    """use group module to add reactions to Kegg pathway

    Args:
        model (libsbml-model): model loaded with libsbml and containing Kegg pathways
        pathway_groups (dict): Kegg Pathway Id as key and reactions Ids as values

    Returns:
        libsbml-model: modified model with groups for pathways
    """
    groups = model.getPlugin('groups')
    group_list = groups.getListOfGroups()
    keys = list(pathway_groups.keys())
    num_reactions = [len(sub) for sub in list(pathway_groups.values())]
    
    for i in range(len(pathway_groups)):
        group_list.createGroup()
        group_list[i].setName(keys[i])
        group_list[i].setId(keys[i]) # important for validity (memote/cobra)
        group_list[i].setMetaId("meta_" + keys[i])
        group_list[i].setKind('partonomy')
        group_list[i].setSBOTerm("SBO:0000633") # NAME
        for j in range(num_reactions[i]):
            group_list[i].createMember()
    
    n_groups = groups.getNumGroups()
    group_list = groups.getListOfGroups()

    for i in range(n_groups):
        group = group_list[i].getName()
        num_members = group_list[i].getNumMembers()
        member_list = group_list[i].getListOfMembers()
        reaction_list = pathway_groups[group]
        for j in range(0, num_members): 
            member_list[j].setIdRef(reaction_list[j])

    return model

def kegg_pathways(modelpath, new_filename):
    """Executes all steps to add KEGG pathways as groups

    Args:
        modelpath (Str): Path to GEM
        new_filename (Str): filename for modified model
    """
    model = load_model_enable_groups(modelpath)

    reactions = extract_kegg_reactions(model)
    pathways = extract_kegg_pathways(reactions)
    pathway_groups = get_pathway_groups(pathways)

    model_pathways = add_kegg_pathways(model, pathways)
    model_pathway_groups = create_pathway_groups(model_pathways, pathway_groups)

    write_to_file(model_pathway_groups, new_filename)