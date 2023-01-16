#!/usr/bin/env python
""" Provides functions for adding KEGG reactions as Group Pathways

If your organism occurs in the KEGG database, extract the KEGG reaction ID from the
annotations of your reactions and identify, in which KEGG pathways this reaction occurs. Add
all KEGG pathways for a reaction then as annotations with the biological qualifier ‘OCCURS_IN’
to the respective reaction.
"""

from tqdm.auto import tqdm
from libsbml import SBMLReader, GroupsExtension
from bioservices import KEGG
from refinegems.load import write_to_file
from refinegems.cvterms import add_cv_term_pathways, parse_id_from_cv_term, add_cv_term_pathways_to_entity

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
    read = reader.readSBMLFromFile(modelpath)  # read from file
    groupextension = GroupsExtension()
    groupURI = groupextension.getURI(3, 1, 1)
    read.enablePackage(groupURI, "groups", True)  # enable groups extension
    read.setPkgRequired('groups', False)  # make groups not required
    model = read.getModel()
    return model


def extract_kegg_reactions(model):
    """extract KEGG Ids from reactions

    Args:
        model (libsbml-model): model loaded with libsbml

    Returns:
        (dict, list): reaction Id as key and Kegg Id as value, Ids of reactions without KEGG annotation
    """
    list_reac = model.getListOfReactions()
    kegg_reactions = {}
    non_kegg_reac = []
    
    for reaction in list_reac:
        kegg_ids = parse_id_from_cv_term(reaction, 'kegg')
        if len(kegg_ids) > 0:
            kegg_reactions[reaction.getId()] = kegg_ids[0]
        else:
            non_kegg_reac.append(reaction.getId())

    return kegg_reactions, non_kegg_reac


def extract_kegg_pathways(kegg_reactions):
    """finds pathway for reactions in model with KEGG Ids
        accesses KEGG API, uses tqdm to report progres to user

    Args:
        kegg_reactions (dict): reaction Id as key and Kegg Id as value

    Returns:
        dict: reaction Id as key and Kegg Pathway Id as value
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
                add_cv_term_pathways_to_entity(path, 'KEGG', reaction)

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


def kegg_pathways(modelpath, new_filename):
    """Executes all steps to add KEGG pathways as groups

    Args:
        modelpath (Str): Path to GEM
        new_filename (Str): filename for modified model
        
    Returns:
        list: Ids of reactions without KEGG annotation
    """
    model = load_model_enable_groups(modelpath)

    reactions, non_kegg_reactions = extract_kegg_reactions(model)
    pathways = extract_kegg_pathways(reactions)
    pathway_groups = get_pathway_groups(pathways)

    model_pathways = add_kegg_pathways(model, pathways)
    model_pathway_groups = create_pathway_groups(
        model_pathways, pathway_groups)

    write_to_file(model_pathway_groups, new_filename)
    
    return non_kegg_reactions
