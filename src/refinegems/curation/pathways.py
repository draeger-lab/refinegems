#!/usr/bin/env python
"""This module provides functions for adding, handling and analysing 
the KEGG pathways (or more specific their annotations) contained in a model.
"""

__author__ = "Famke Baeuerle and Carolin Brune"

############################################################################
# requirements
############################################################################

import cobra
import re

from tqdm.auto import tqdm
from libsbml import SBMLReader, GroupsExtension
from libsbml import Model as libModel
from bioservices import KEGG
from ..utility.cvterms import add_cv_term_pathways, get_id_from_cv_term, add_cv_term_pathways_to_entity
from ..classes.reports import KEGGPathwayAnalysisReport

############################################################################
# functions
############################################################################

# adding KEGG reactions as Group Pathways
# ---------------------------------------
""" 
If your organism occurs in the KEGG database, 
extract the KEGG reaction ID from the annotations of your reactions and identify, 
in which KEGG pathways this reaction occurs. 
Add all KEGG pathways for a reaction then as annotations 
with the biological qualifier `OCCURS_IN` to the respective reaction.
"""

def load_model_enable_groups(modelpath: str) -> libModel:
    """Loads model as document using libSBML and enables groups extension

    Args:
        - modelpath (str): 
            Path to GEM

    Returns:
        libModel: 
            Model loaded with libSBML
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    groupextension = GroupsExtension()
    groupURI = groupextension.getURI(3, 1, 1)
    read.enablePackage(groupURI, "groups", True)  # enable groups extension
    read.setPkgRequired('groups', False)  # make groups not required
    model = read.getModel()
    return model


def extract_kegg_reactions(model: libModel) -> tuple[dict, list]:
    """Extract KEGG Ids from reactions

    Args:
        - model (libModel): 
            Model loaded with libSBML. Output of :py:func:`~refinegems.curation.pathways.load_model_enable_groups`.

    Returns:
        tuple: 
            Dictionary 'reaction_id': 'KEGG_id' (1) & List of reactions without KEGG Id (2)
            
            (1) dict: Reaction Id as key and KEGG Id as value
            (2) list: Ids of reactions without KEGG annotation
    """
    list_reac = model.getListOfReactions()
    kegg_reactions = {}
    non_kegg_reac = []
    
    for reaction in list_reac:
        kegg_ids = get_id_from_cv_term(reaction, 'kegg.reaction')
        if len(kegg_ids) > 0:
            kegg_reactions[reaction.getId()] = kegg_ids[0]
        else:
            non_kegg_reac.append(reaction.getId())

    return kegg_reactions, non_kegg_reac


def extract_kegg_pathways(kegg_reactions: dict) -> dict:
    """Finds pathway for reactions in model with KEGG Ids, accesses KEGG API, uses tqdm to report progres to user

    Args:
        - kegg_reactions (dict): 
            Reaction Id as key and KEGG Id as value. Output[0] from :py:func:`~refinegems.curation.pathways.extract_kegg_reactions`.

    Returns:
        dict: 
            Reaction Id as key and KEGG Pathway Id as value
    """
    k = KEGG()
    kegg_pathways = {}

    print('Extracting pathway Id for each reaction:')
    for reaction in tqdm(kegg_reactions.keys()):
        kegg_reaction = k.get(kegg_reactions[reaction])
        dbentry = k.parse(kegg_reaction)
        # sometimes parse does not work -> try and except
        try:
            pathways = [x for x in dbentry['PATHWAY']]
        except BaseException:
            pathways = []
        kegg_pathways[reaction] = pathways

    return kegg_pathways


def add_kegg_pathways(model, kegg_pathways) -> libModel:
    """Add KEGG reactions as BQB_OCCURS_IN.

    Args:
        - model (libModel): 
            Model loaded with libSBML. Output of :py:func:`~refinegems.curation.pathways.load_model_enable_groups`.
        - kegg_pathways (dict): 
            Reaction Id as key and KEGG Pathway Id as value. Output of :py:func:`~refinegems.curation.pathways.extract_kegg_pathways`.

    Returns:
        libModel: 
        Modified model with KEGG pathways.
    """
    list_reac = model.getListOfReactions()

    for reaction in list_reac:
        if reaction.getId() in kegg_pathways.keys():
            for path in kegg_pathways[reaction.getId()]:
                add_cv_term_pathways_to_entity(path, 'KEGG', reaction)

    return model


def get_pathway_groups(kegg_pathways) -> dict:
    """Group reaction into pathways.

    Args:
        - kegg_pathways (dict): 
            Reaction Id as key and KEGG Pathway Id as value. Output of :py:func:`~refinegems.curation.pathways.extract_kegg_pathways`.

    Returns:
        dict: 
            KEGG Pathway Id as key and reactions Ids as values.
    """
    pathway_groups = {}
    for reaction in kegg_pathways.keys():
        for path in kegg_pathways[reaction]:
            if path not in pathway_groups.keys():
                pathway_groups[path] = [reaction]
            else:
                pathway_groups[path].append(reaction)
    return pathway_groups


def create_pathway_groups(model: libModel, pathway_groups) -> libModel:
    """Use group module to add reactions to KEGG pathway.

    Args:
        - model (libModel): 
            Model loaded with libSBML. Output of :py:func:`~refinegems.curation.pathways.load_model_enable_groups`.
        - pathway_groups (dict): 
            KEGG Pathway Id as key and reactions Ids as values. Output of :py:func:`~refinegmes.curation.pathways.get_pathway_groups`.

    Returns:
        libModel: 
            Modified model with groups for pathways.
    """
    k = KEGG()
    groups = model.getPlugin('groups')
    group_list = groups.getListOfGroups()
    keys = list(pathway_groups.keys())
    num_reactions = [len(sub) for sub in list(pathway_groups.values())]

    print('Adding pathways as groups to the model:')
    for i in tqdm(range(len(pathway_groups))):
        kegg_pathway = k.get(keys[i])
        dbentry = k.parse(kegg_pathway)
        if groups.getGroup('G_' + keys[i]) is not None:
            group = groups.getGroup('G_' + keys[i])
            group.setName(dbentry['NAME'][0])
            group.setMetaId("meta_" + 'G_' + keys[i])
            group.setKind('partonomy')
            group.setSBOTerm("SBO:0000633")  # NAME
            add_cv_term_pathways(keys[i], 'KEGG', group)
            for reac in pathway_groups[keys[i]]:
                if group.getMemberByIdRef(reac) is None:
                    member = group.createMember()
                    member.setIdRef(reac)
        else:
            group = group_list.createGroup() 
            group.setName(dbentry['NAME'][0])
            group.setId('G_' + keys[i])  # important for validity (memote/cobra)
            group.setMetaId("meta_" + 'G_' + keys[i])
            group.setKind('partonomy')
            group.setSBOTerm("SBO:0000633")  # NAME
            add_cv_term_pathways(keys[i], 'KEGG', group)
            for reac in pathway_groups[keys[i]]:
                if group.getMemberByIdRef(reac) is None:
                    member = group.createMember()
                    member.setIdRef(reac)

    return model


def kegg_pathways(modelpath: str) -> tuple[libModel, list[str]]:
    """Executes all steps to add KEGG pathways as groups

    Args:
        - modelpath (str): 
            Path to GEM.
        
    Returns:
        tuple: 
            libSBML model (1) & List of reactions without KEGG Id (2)

            (1) libModel: Modified model with Pathways as groups
            (2) list: Ids of reactions without KEGG annotation
    """
    model = load_model_enable_groups(modelpath)

    reactions, non_kegg_reactions = extract_kegg_reactions(model)
    pathways = extract_kegg_pathways(reactions)
    pathway_groups = get_pathway_groups(pathways)

    model_pathways = add_kegg_pathways(model, pathways)
    model_pathway_groups = create_pathway_groups(
        model_pathways, pathway_groups)
    
    return model_pathway_groups, non_kegg_reactions


# analyse the pathways in a model
# -------------------------------

def kegg_pathway_analysis(model:cobra.Model) -> KEGGPathwayAnalysisReport:
    """Analyse the pathways that are covered by the model.

    The analysis is based on the KEGG pathway classification and the available
    KEGG pathway identifiers present in the model.
    
    Note: one reaction can have multiple pathway identifiers associated with it. 
    This analysis focuses on the total number of IDs found within the model.

    Args:
        - model (cobra.Model): 
            A model loaded with COBRApy.

    Returns:
        KEGGPathwayAnalysisReport: 
            The KEGG pathway analysis report.
    """
    # create report
    report = KEGGPathwayAnalysisReport(total_reac=len(model.reactions))

    pathways = dict()
    counter = 0
    # extract KEGG pathway IDs from all reactions
    for r in model.reactions:
        if 'kegg.pathway' in r.annotation.keys():
            counter += 1
            anno = r.annotation['kegg.pathway']
            # case 1: only one annotation found
            if isinstance(anno,str):
                anno = re.sub(r'^[a-z]*','',anno)
                if anno in pathways:
                    pathways[anno] += 1
                else:
                    pathways[anno] = 1
            # case 2: multiple annotations for one reaction found
            else:
                for x in anno:
                    x = re.sub(r'^[a-z]*','',x)
                    if x in pathways:
                        pathways[x] += 1
                    else:
                        pathways[x] = 1
    
    # add counter to report
    report.kegg_count = counter

    # identify global and overview pathway identifier
    global_map = {}
    over_map = {}
    rest = {}
    for k,v in pathways.items():
        if k.startswith('011'):
            global_map[k] = v
        elif k.startswith('012'):
            over_map[k] = v
        else:
            rest[k] = v

    # add IDs in corresponding class to report
    report.kegg_global = global_map
    report.kegg_over = over_map
    report.kegg_paths = rest

    return report
