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
import urllib

from bioservices import KEGG
from Bio.KEGG import REST, Enzyme
from libsbml import SBMLReader, GroupsExtension
from libsbml import Model as libModel
from tqdm.auto import tqdm

from ..utility.cvterms import (
    add_cv_term_pathways,
    get_id_from_cv_term,
    add_cv_term_pathways_to_entity,
)
from ..utility.db_access import kegg_reaction_parser
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
    read.setPkgRequired("groups", False)  # make groups not required
    model = read.getModel()
    return model


def _extract_kegg_ec_from_reac(model: libModel) -> tuple[dict, list]:
    """Given a model. extract KEGG IDs and EC numbers from all reactions.

    Args:
        - model (libModel):
            A model loaded with libSBML.

    Returns:
        tuple[dict, list]:
            Tuple with (1) a dictionary with the mappings and (2) a list of unmapped reactions.

            (1) Dictionary:
                Reaction Id as key and dictionary as value.
                The dictionary contains either one or both of the following
                keys with the corresponding values:

                - 'kegg.reaction': List of KEGG reaction Ids
                - 'eccode': List of EC numbers

            (2) List:
                Reaction Ids without KEGG or EC annotation.
    """

    list_reac = model.getListOfReactions()
    mapped = {}
    non_mapped = []

    for reaction in list_reac:
        ecs = get_id_from_cv_term(reaction, "ec-code")
        kegg_ids = get_id_from_cv_term(reaction, "kegg.reaction")
        mappings = dict()
        if len(kegg_ids) > 0:
            mappings["kegg.reaction"] = kegg_ids
        if len(ecs) > 0:
            mappings["eccode"] = ecs
        if len(mappings) == 0:
            non_mapped.append(reaction.getId())
        else:
            mapped[reaction.getId()] = mappings

    # mapped : dict = {'reac_id': {kegg: [], ec-code:[]}}
    return mapped, non_mapped


def find_kegg_pathways(
    mapped_reacs: dict, viaEC: bool = False, viaRC: bool = False
) -> dict:
    """Given a dictionary of reaction IDs mapped to KEGG reaction IDs and/or EC numbers,
    extract the KEGG pathways for each reaction based on the KEGG reaction ID.

    Args:
        - mapped_reacs (dict):
            Dictionary containing the information about the reactions.
            For more information see, :py:func:`~refinegems.curation.pathways._extract_kegg_ec_from_reac`.
        - viaEC (bool, optional):
            If True, also tries mapping to pathways via EC number, if
            via reaction ID is unsuccessful.
            Defaults to False.
        - viaRC (bool, optional):
            If True, also tries mapping to pathways via reaction class, if
            via reaction ID is unsuccessful.
            Defaults to False.

    Returns:
        dict:
            Dictionary with the reaction IDs as keys and a list of KEGG pathway IDs as values.
    """

    def _get_pathway_via_rc(rc_list: list[str]) -> list[str]:
        """Helper function to extract KEGG pathways IDs by a list of reaction class IDs.

        Args:
            - rc_list (list[str]):
                List of reaction class IDs.

        Returns:
            list[str]:
                List of extracted pathways.
                If none are found, return an empty list.
        """

        pathway_ids = []
        collect = False
        for kegg_rc in rc_list:
            kegg_rc = REST.kegg_get(kegg_rc)
            kegg_rc = kegg_rc.read()
            for line in kegg_rc.split("\n"):
                if line:
                    if line.startswith("PATHWAY"):
                        collect = True
                        pathway_ids.append(
                            line.replace("PATHWAY", "", 1).strip().split(" ")[0]
                        )
                    elif collect == True and line[0] != "/":
                        if line[0].isupper():
                            collect = False
                        else:
                            pathway_ids.append(line.strip().split(" ")[0])

        return pathway_ids

    kegg_pathways = {}

    print("Extracting pathway Id for each reaction:")
    for reac_id in tqdm(mapped_reacs.keys(), unit="reaction"):

        kegg_reaction = None
        pathways = []

        # via reaction
        if "kegg.reaction" in mapped_reacs[reac_id].keys():
            for kegg_id in mapped_reacs[reac_id]["kegg.reaction"]:
                try:
                    kegg_reaction = kegg_reaction_parser(kegg_id)
                    if isinstance(kegg_reaction["db"]["kegg.pathway"], list):
                        pathways.extend(kegg_reaction["db"]["kegg.pathway"])
                    else:
                        pathways.append(kegg_reaction["db"]["kegg.pathway"])

                # exception handling
                except urllib.error.HTTPError:
                    print(f"HTTPError: {reac_id} on {kegg_id}")
                except ConnectionResetError:
                    print(f"ConnectionResetError: {reac_id} on {kegg_id}")
                except urllib.error.URLError:
                    print(f"URLError: {reac_id} on {kegg_id}")
                except KeyError:
                    # no pathway found
                    pass
                # except Exception as e:
                #     print(F'Something unexpected happened: {reac_id} on {kegg_id}')
                #     print(repr(e))

        # via RC
        # via reaction class
        # can lead to some additional classes, as the RC are not as strictly defined as
        # the reactions themselves

        if viaRC and kegg_reaction and len(pathways) == 0:

            try:
                rc_results = _get_pathway_via_rc(kegg_reaction["db"]["kegg.rclass"])
                if len(rc_results) > 0:
                    pathways.extend(rc_results)
            # exception handling
            except urllib.error.HTTPError:
                print(f"HTTPError: {reac_id} on RCLASS")
            except ConnectionResetError:
                print(f"ConnectionResetError: {reac_id} on RCLASS")
            except urllib.error.URLError:
                print(f"URLError: {reac_id} on RCLASS")
            except KeyError:
                # no pathway found
                pass
            except Exception as e:
                print(f"Something unexpected happened: {reac_id} on RCLASS")

        # via EC
        # seems really sketchy to do it this way, as ONE EC number
        # can include MANY reactions for different pathways

        if viaEC and len(pathways) == 0:
            ec = None
            # tests, if reaction already has an annotated EC number
            if "eccode" in mapped_reacs[reac_id].keys():
                ec = mapped_reacs[reac_id]["eccode"]
            # try to get via kegg reaction
            elif (
                kegg_reaction
                and "db" in kegg_reaction.keys()
                and "ec-code" in kegg_reaction["db"]
            ):
                ec = kegg_reaction["db"]["ec-code"]

            # if EC number found,
            # get pathways information using the EC number
            if ec:
                ecnum = None
                try:
                    for ecnum in ec:
                        kegg_ec = REST.kegg_get(f"ec:{ecnum}")
                        kegg_ec = Enzyme.read(kegg_ec)
                        # no pathway found
                        if len(kegg_ec.pathway) == 0 or kegg_ec.pathway == None:
                            pass
                        # pathway found
                        else:
                            for i in kegg_ec.pathway:
                                pathways.append(i[1])
                # exception handling
                except urllib.error.HTTPError:
                    print(f"HTTPError: {reac_id} on {ecnum}")
                except ConnectionResetError:
                    print(f"ConnectionResetError: {reac_id} on {ecnum}")
                except urllib.error.URLError:
                    print(f"URLError: {reac_id} on {ecnum}")
                except KeyError:
                    # no pathway found
                    pass
                except Exception as e:
                    print(f"Something unexpected happened: {reac_id} on {ecnum}")

        # store the pathways
        kegg_pathways[reac_id] = list(set(pathways))

    return kegg_pathways


def add_kegg_pathways(model, kegg_pathways) -> libModel:
    """Add KEGG reactions as BQB_OCCURS_IN.

    Args:
        - model (libModel):
            Model loaded with libSBML. Output of :py:func:`~refinegems.curation.pathways.load_model_enable_groups`.
        - kegg_pathways (dict):
            Reaction Id as key and KEGG Pathway Id as value, e.g. see output of :py:func:`~refinegems.curation.pathways.find_kegg_pathways`.

    Returns:
        libModel:
        Modified model with KEGG pathways.
    """
    list_reac = model.getListOfReactions()

    for reaction in list_reac:
        if reaction.getId() in kegg_pathways.keys():
            for path in kegg_pathways[reaction.getId()]:
                add_cv_term_pathways_to_entity(path, "KEGG", reaction)

    return model


def create_pathway_groups(model: libModel, pathway_groups) -> libModel:
    """Use group module to add reactions to KEGG pathway.

    Args:
        - model (libModel):
            Model loaded with libSBML. Output of :py:func:`~refinegems.curation.pathways.load_model_enable_groups`.
        - pathway_groups (dict):
            KEGG Pathway Id as key and reactions Ids as values, e.g. see output of :py:func:`~refinegmes.curation.pathways.set_kegg_pathways._invert_reac_pathway_dict`.

    Returns:
        libModel:
            Modified model with groups for pathways.
    """
    k = KEGG()
    groups = model.getPlugin("groups")
    group_list = groups.getListOfGroups()
    keys = list(pathway_groups.keys())

    print("Adding pathways as groups to the model:")
    for i in tqdm(range(len(pathway_groups))):
        kegg_pathway = k.get(keys[i])
        dbentry = k.parse(kegg_pathway)
        if groups.getGroup("G_" + keys[i]) is not None:
            group = groups.getGroup("G_" + keys[i])
            group.setName(dbentry["NAME"][0])
            group.setMetaId("meta_" + "G_" + keys[i])
            group.setKind("partonomy")
            group.setSBOTerm("SBO:0000633")  # NAME
            add_cv_term_pathways(keys[i], "KEGG", group)
            for reac in pathway_groups[keys[i]]:
                if group.getMemberByIdRef(reac) is None:
                    member = group.createMember()
                    member.setIdRef(reac)
        else:
            group = group_list.createGroup()
            group.setName(dbentry["NAME"][0])
            group.setId("G_" + keys[i])  # important for validity (memote/cobra)
            group.setMetaId("meta_" + "G_" + keys[i])
            group.setKind("partonomy")
            group.setSBOTerm("SBO:0000633")  # NAME
            add_cv_term_pathways(keys[i], "KEGG", group)
            for reac in pathway_groups[keys[i]]:
                if group.getMemberByIdRef(reac) is None:
                    member = group.createMember()
                    member.setIdRef(reac)

    return model


def set_kegg_pathways(
    modelpath: str, viaEC: bool = False, viaRC: bool = False
) -> tuple[libModel, list[str]]:
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

    def _invert_reac_patway_dict(kegg_pathways) -> dict:
        """Group reaction into pathways.

        Args:
            - kegg_pathways (dict):
                Reaction Id as key and KEGG Pathway Id as value. Output of :py:func:`~refinegems.curation.pathways.find_kegg_pathways`.

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

    # load model with groups enabled
    model = load_model_enable_groups(modelpath)

    # extract information about KEGG and EC numbers from model reactions
    reactions, non_kegg_reactions = _extract_kegg_ec_from_reac(model)

    # @DEBUG .................
    # reactions = {k:reactions[k] for idx,k in enumerate(reactions) if idx < 5}
    # ........................

    # add kegg pathways
    pathways = find_kegg_pathways(reactions, viaEC=viaEC, viaRC=viaRC)
    model_pathways = add_kegg_pathways(model, pathways)

    # add corresponding groups
    pathway_groups = _invert_reac_patway_dict(pathways)
    model_pathway_groups = create_pathway_groups(model_pathways, pathway_groups)

    return model_pathway_groups, non_kegg_reactions


# analyse the pathways in a model
# -------------------------------


def kegg_pathway_analysis(model: cobra.Model) -> KEGGPathwayAnalysisReport:
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
        if "kegg.pathway" in r.annotation.keys():
            counter += 1
            anno = r.annotation["kegg.pathway"]
            # case 1: only one annotation found
            if isinstance(anno, str):
                anno = re.sub(r"^[a-z]*", "", anno)
                if anno in pathways:
                    pathways[anno] += 1
                else:
                    pathways[anno] = 1
            # case 2: multiple annotations for one reaction found
            else:
                for x in anno:
                    x = re.sub(r"^[a-z]*", "", x)
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
    for k, v in pathways.items():
        if k.startswith("011"):
            global_map[k] = v
        elif k.startswith("012"):
            over_map[k] = v
        else:
            rest[k] = v

    # add IDs in corresponding class to report
    report.kegg_global = global_map
    report.kegg_over = over_map
    report.kegg_paths = rest

    return report
