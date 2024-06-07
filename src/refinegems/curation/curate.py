#!/usr/bin/env python
"""General functions for curating a model.

Includes functions for adding annotations from manual tables, annotations synchronisation, 
handling duplicates and more."""

__author__ = "Famke Baeuerle und Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import pandas as pd
import re
import warnings

from libsbml import Model as libModel
from tqdm.auto import tqdm
from typing import Literal

from ..utility.cvterms import add_cv_term_reactions, add_cv_term_metabolites, metabol_db_dict, get_id_from_cv_term
from ..utility.entities import create_gpr_from_locus_tag, create_reaction
from .biomass import test_biomass_presence

################################################################################
# variables
################################################################################

NH_PATTERN = re.compile('nh[3-4]')

################################################################################
# functions
################################################################################
    
# add information from manual curation
# ------------------------------------

# @ TODO : check against gapfill & co for new-ish functionalities and delete        
def add_reactions_from_table(model: libModel, table: pd.DataFrame, email: str) -> libModel:
    """Wrapper function to use with table format given in data/manual_curation.xlsx, sheet gapfill: Adds all reactions with their info given in the table to the given model

    **Deprecation warning**: will be deprecated in a future update.

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - table (pd-DataFrame): 
            Table in format of sheet gapfill from manual_curation.xlsx located in the data folder
        - email (str): 
            User Email to access the NCBI Entrez database

    Returns:
        libModel: 
            Modified model with new reactions
    """
    mes = f'add_reactions_from_table will be deprecated in a future update.'
    warnings.warn(mes,type=FutureWarning)

    for reaction_info in tqdm(table.groupby('BIGG')):
        reac_id = reaction_info[0]
        if model.getReaction(str(reac_id)) is None:
            reactants = dict(table.loc[table['BIGG'] == reac_id, ['educts', 'stoich_e']].dropna().values)
            products = dict(table.loc[table['BIGG'] == reac_id, ['products', 'stoich_p']].dropna().values)
            fluxes = table.loc[table['BIGG'] == reac_id, ['lower_bound', 'upper_bound']].dropna().to_dict('records')[0]
            name = table.loc[table['BIGG'] == reac_id, ['name']].dropna().to_dict('records')[0]['name']
            reversible = table.loc[table['BIGG'] == reac_id, ['reversible']].dropna().to_dict('records')[0]['reversible']
            fast = table.loc[table['BIGG'] == reac_id, ['fast']].dropna().to_dict('records')[0]['fast']
            try:
                sbo = table.loc[table['BIGG'] == reac_id, ['sbo']].dropna().to_dict('records')[0]['sbo']
            except (IndexError):
                print('SBO Term for ' + str(reac_id) + ' will be set to standard "SBO:0000167" (biochemical or transport reaction)')
                sbo = "SBO:0000167"
            reaction, model = create_reaction(model, reac_id, name, reactants, products, fluxes, reversible, fast, sbo)
            for (columnName, columnData) in table.loc[table['BIGG'] == reac_id].drop(['educts', 'stoich_e','products', 'stoich_p','lower_bound', 'upper_bound', 'name', 'reversible', 'fast'], axis=1).fillna(0).iteritems():
                for entry in columnData.values:
                    if not entry == 0:
                        if columnName == 'locus':
                            reaction.getPlugin(0).createGeneProductAssociation().createGeneProductRef().setGeneProduct(str(entry))
                            if model.getPlugin(0).getGeneProductByLabel(str(entry)) is None:
                                gpr, model = create_gpr_from_locus_tag(model, str(entry), email)
                        else:
                            add_cv_term_reactions(str(entry), str(columnName), reaction)
    return model


# @ TODO : check against gapfill & co for new-ish functionalities and delete     
def update_annotations_from_table(model: libModel, table: pd.DataFrame) -> libModel:
    """Wrapper function to use with table format given in data/manual_curation.xlsx, sheet metabs: Updates annotation of metabolites given in the table

    **Deprecation warning**: will be deprecated in a future update.

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - table (pd-DataFrame): 
            Table in format of sheet metabs from manual_curation.xlsx located in the data folder

    Returns:
        libModel: 
            Modified model with new annotations
    """
    mes = f'update_annotations_from_table will be deprecated in a future update.'
    warnings.warn(mes,type=FutureWarning)

    table = table.drop(['Name', 'FORMULA', 'Notiz'], axis=1).fillna(0)
    table['PUBCHEM'] = table['PUBCHEM'].astype(int)
    for metab_info in tqdm(table.groupby('BIGG')):
        met = metab_info[0]
        for comp in ['_c', '_e', '_p']:
            try:
                metab = model.getSpecies('M_' + met + comp)
                #metab.unsetAnnotation()
                if not metab.isSetMetaId():
                    metab.setMetaId('meta_' + metab.getId())
                for (columnName, columnData) in table.loc[table['BIGG'] == met].iteritems():
                    for entry in columnData.values:
                        if not entry == 0:
                            add_cv_term_metabolites(str(entry), str(columnName), metab)
            except (AttributeError):
                print(met + comp + ' not in model')
    return model


# sychronise annotations
# ----------------------

def update_annotations_from_others(model: libModel) -> libModel:
    """Synchronizes metabolite annotations for core, periplasm and extracelullar

    Args:
        - model (libModel): 
            Model loaded with libSBML

    Returns:
        libModel: 
            Modified model with synchronized annotations
    """
    for metab in model.getListOfSpecies():
        base = metab.getId()[:-2]
        for comp in ['_c', '_e', '_p']:
            other_metab = model.getSpecies(base + comp)
            if other_metab is not None:
                if not other_metab.isSetMetaId():
                    other_metab.setMetaId('meta_' + other_metab.getId())
                for db_id, code in metabol_db_dict.items():
                    id = get_id_from_cv_term(metab, code)
                    for entry in id:
                        if entry is not None:
                            add_cv_term_metabolites(entry, db_id, other_metab)    
    return model               


def complete_BioMetaCyc(model:cobra.Model) -> cobra.Model:
    """Check for existing MetaCyc and BioCyc annotations for metabolites and
    reactions and generate them from the other if one of the two is missing.

    Args:
        model (cobra.Model): 
            A genome-scale model to be checked 
            for complete BioCyc/MetaCyc annotations.

    Returns:
        cobra.Model: 
            The updated model.
    """


    # reactions
    # ---------
    for reac in model.reactions:
        # case 1: MetaCyc but not BioCyc
        if 'metacyc.reaction' in reac.annotation and not 'biocyc' in reac.annotation:
            # add database (META) information to get BioCyc ID
            if isinstance(reac.annotation['metacyc.reaction'], list):
                reac.annotation['biocyc'] = ['META:' + _ for _ in reac.annotation['metacyc.reaction']]
            else:
                reac.annotation['biocyc'] = 'META:' + reac.annotation['metacyc.reaction']
        # case 2: BioCyc, but no MetaCyc
        elif 'biocyc' in reac.annotation and not 'metacyc.reaction' in reac.annotation:
            # if there are multiple metacyc.reaction annotation
            if isinstance(reac.annotation['biocyc'], list):
                add_anno = []
                for biocyc_anno in reac.annotation['biocyc']:
                    if ':' in biocyc_anno:
                        # exclude organism identifier from ID to get MetaCyc ID
                        add_anno.append(biocyc_anno.split(':')[1])
                    else:
                        # high possibility that information is faulty - do not use it
                        print(F'\n\nWarning: Unexpected BioCyc annotation {biocyc_anno} for reaction {reac.id}')
                reac.annotation['metacyc.reaction'] = add_anno
            # if there is only one
            else:
                if ':' in reac.annotation['biocyc']:
                    # exclude organism identifier from ID to get MetaCyc ID
                    reac.annotation['metacyc.reaction'] = reac.annotation['biocyc'].split(':')[1]
                else:
                    # high possibility that information is faulty - do not use it
                    print(F'\n\nWarning: Unexpected BioCyc annotation {reac.annotation["biocyc"]} for reaction {reac.id}')
        # case 3: both or no = skip
        else:
            continue

    # metabolites
    # -----------
    for meta in model.metabolites:
        # case 1: MetaCyc but not BioCyc
        if 'metacyc.compound' in meta.annotation and not 'biocyc' in meta.annotation:
            # add database (META) information to get BioCyc ID
            if isinstance(meta.annotation['metacyc.compound'],list):
                meta.annotation['biocyc'] = ['META:' + _ for _ in meta.annotation['metacyc.compound']]
            else:
                meta.annotation['biocyc'] = 'META:' + meta.annotation['metacyc.compound']
        # case 2: BioCyc, but no MetaCyc
        elif 'biocyc' in meta.annotation and not 'metacyc.compound' in meta.annotation:
            # if there are multiple metacyc.compound annotations
            if isinstance(meta.annotation['biocyc'], list):
                add_anno = []
                for biocyc_anno in meta.annotation['biocyc']:
                    if ':' in biocyc_anno:
                        # exclude organism identifier from ID to get MetaCyc ID
                        add_anno.append(biocyc_anno.split(':')[1])
                    else:
                        # high possibility that information is faulty - do not use it
                        print(F'\n\nWarning: Unexpected BioCyc annotation {biocyc_anno} for metabolite {meta.id}')
                meta.annotation['metacyc.compound'] = add_anno
            # if there is only one
            else:
                if ':' in meta.annotation['biocyc']:
                    # exclude organism identifier from ID to get MetaCyc ID
                    meta.annotation['metacyc.compound'] = meta.annotation['biocyc'].split(':')[1]
                else:
                    # high possibility that information is faulty - do not use it
                    print(F'\n\nWarning: Unexpected BioCyc annotation {meta.annotation["biocyc"]} for metabolite {meta.id}')
        # case 3: both or no = skip
        else:
            continue

    return model



# duplicates
# ----------

def resolve_duplicate_reactions(model:cobra.Model, based_on:str='reaction', remove_reac:bool=True) -> cobra.Model:
    """Resolve and remove duplicate reaction based on their reaction equation
    and matching database identifiers. Only if all match or are nan will one of
    the reactions be removed.

    Args:
        - model (cobra.Model): 
            A model loaded with COBRApy.
        - based_on (str, optional): 
            Label to base the resolvement process on . 
            Can be 'reaction' or any other annotation label. 
            Defaults to 'reaction'.
        - remove_reac (bool, optional): 
            When True, combines and remove duplicates. 
            Otherwise only reports the findings.
            Defaults to True.

    Returns:
        cobra.Model: 
            The model.
    """

    # get annotation and compartment information
    anno_reac = []
    for r in model.reactions:
        anno_reac.append({'id':r.id, 'compartment':str(r.compartments), 'reaction':r.reaction} | r.annotation)
    df_reac = pd.DataFrame.from_dict(anno_reac)

    # check if based_on is valid
    if not based_on in df_reac.columns.tolist():
        warnings.warn(F'Warning: Annotation column {based_on} does not exists. Search for duplicates will be skipped.')
        return model

    # set basic parameters
    skip_cols = ['id','compartment','bigg.reaction','reaction',based_on]
    colnames = df_reac.columns.tolist()

    for c in df_reac.groupby('compartment'):
        # note: using groupby drops nans
        for mnx in c[1].groupby(based_on):
            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                # check annotations
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = mnx[1][col].dropna().explode().unique().tolist()
                    else:
                        dupl = False
                        break

                # if duplicate found
                if dupl:
                    if remove_reac:

                        # choose reaction to keep
                        keep_reac = model.reactions.get_by_id(mnx[1]['id'].tolist()[0])

                        # resolve annotations
                        for key, value in annotations.items():
                            if len(value) > 0 and not key in keep_reac.annotation:
                                keep_reac.annotation[key] = value

                        # combine gene reaction rules
                        for r_id in mnx[1]['id'].tolist()[1:]:
                            keep_reac.gene_reaction_rule = keep_reac.gene_reaction_rule + model.reactions.get_by_id(r_id).gene_reaction_rule
                            model.reactions.get_by_id(r_id).delete()
                            print(F'\tDuplicate reaction {r_id} found. Combined to {keep_reac.id} and deleted.')
                    else:
                        print(f'\tDuplicate reactions {", ".join(mnx[1]["id"].tolist())} found.')

    return model


def resolve_duplicate_metabolites(model:cobra.Model, based_on:str='metanetx.chemical', replace:bool=True) -> cobra.Model:
    """Resolve duplicate metabolites in a model. Metabolites are considered
    duplicate if they share the same annotations (same or nan).
    Note: Depending on the starting database, the results might differ.

    Args:
        - model (cobra.Model): 
            The model loaded with COBRApy.
        - based_on (str, optional): 
            Label to base the resolvement process on . 
            Can be any annotation label. 
            Defaults to 'metanetx.chemical'.
        - replace (bool, optional): 
            Either report the duplicates (False) 
            or replace them with one (True). Defaults to True.

    Returns:
        cobra.Model: 
            The model.
    """

    # get annotation and compartment information
    anno_meta = []
    for m in model.metabolites:
        anno_meta.append({'id':m.id, 'compartment':m.compartment} | m.annotation)
    df_meta = pd.DataFrame.from_dict(anno_meta)

    # check if based_on is valid
    if not based_on in df_meta.columns.tolist():
        warnings.warn(F'Warning: Annotation {based_on} not found. Search for metabolite duplicates skipped.')
        return model

    # set basic parameters
    skip_cols = ['id','compartment','bigg.metabolite',based_on]
    colnames = df_meta.columns.tolist()
    # get objective function
    bof_list = test_biomass_presence(model)
    if len(bof_list) == 1:
        objective_function = bof_list[0]
    elif len(bof_list) > 1:
        mes = f'Multiple BOFs detected. Will be using {bof_list[0]}'
        warnings.warn(mes,category=UserWarning)
        objective_function = bof_list[0]
    else:
        mes = 'No BOF detected. Might lead to problems during duplicate removal.'
        warnings.warn(mes, category=UserWarning) 

    for c in df_meta.groupby('compartment'):
        # note: using groupby drops nans
        # note: mnx as starting point was choosen, as it seems to have the most annotations (easy to get)
        for mnx in c[1].groupby(based_on):

            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    # check annotation, if current group truly consists of duplicates
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = mnx[1][col].dropna().explode().unique().tolist()
                    else:
                        dupl = False
                        break
            # no duplicates = check next
            else:
                continue

            # if duplicate is found:
            if dupl:
                # either replace ...
                if replace:
                    # choose metabolite to keep
                    keep_meta = model.metabolites.get_by_id(mnx[1]['id'].tolist()[0])
                    # resolve annotations
                    for key, value in annotations.items():
                        if len(value) > 0 and not key in keep_meta.annotation:
                            keep_meta.annotations[key] = value
                    # note: charge and formula should be in valid range and be corrected by MCC (if needed)

                    # replace duplicates with the metabolite to be kept
                    #     to ensure consistency, only delete duplicate metabolites, which
                    #     do NOT share ANY reactions

                    # retrieve reactions for metabolites set for keeping
                    keep_reac = [_.id for _ in model.metabolites.get_by_id(keep_meta.id).reactions]
                    # iterate over metabolites set for deletion
                    for del_meta_id in mnx[1]['id'].tolist()[1:]:
                        # retrieve reaction for metabolites set for deletion
                        del_reac = [_.id for _ in model.metabolites.get_by_id(del_meta_id).reactions]
                        # get intersection of reactions (keep + del)
                        reac_intersec = list(set(keep_reac) & set(del_reac))

                        # if intersection empty, metabolite is with a high probability indeed a duplicate
                        # Special case: NH3 / NH4
                        #    intersection does not have to be emtpy, know 'problem' caused by CarveMe
                        if len(reac_intersec) == 0 or all([re.search(NH_PATTERN,_) for _ in [keep_meta.id,del_meta_id]]):
                            # automated deletion is only advisable, if consistency can be retained
                            perform_deletion = True
                            with model as model_del:

                                # if the special case is detected ...
                                if all([re.search(NH_PATTERN,_) for _ in [keep_meta.id,del_meta_id]]):
                                    print(F'\tSpecial case -Duplicate NH4/NH3- detected.\n\tTrying to solve by additionally removing reactions containing both metabolites.')
                                    # ... remove reactions with nh3 and nh4 both present
                                    for del_reac_id in reac_intersec:
                                        # if objective_function is part of the set
                                        # automated deletion is (currently) not possible
                                        if del_reac_id == objective_function:
                                            perform_deletion = False
                                            break
                                        model_del.reactions.get_by_id(del_reac_id).remove_from_model()
                                    # set the metabolites to be deleted to be the one NOT in the objective functions
                                    # to avoid inconsistencies
                                    if del_meta_id in [_.id for _ in model_del.reactions.get_by_id(objective_function).metabolites]:
                                        temp = del_meta_id
                                        del_meta_id = keep_meta.id
                                        keep_meta = model_del.metabolites.get_by_id(temp)

                                # try replacing metabolite with the kept duplicate ...
                                for reac in model_del.metabolites.get_by_id(del_meta_id).reactions:

                                    reac.reaction = reac.reaction.replace(del_meta_id,keep_meta.id)

                                    # skip if objective function is found
                                    if reac.id == objective_function:
                                        continue

                                    # check if consistency is still intact
                                    balance_test = reac.check_mass_balance()
                                    if not reac.boundary and len(balance_test) > 0:
                                        # try fixing H-balance
                                        if 'H' in balance_test.keys():
                                            # ..............................
                                            # TODO:
                                            #    get H according to compartment
                                            #    current implementation relies heavily
                                            #    on 'correct' use input: compartment should have format C_c or c (C_p, p, C_e, e etc.)
                                            # ..............................
                                            reac_comp = reac.compartments.pop()[-1]
                                            if reac_comp == 'c':
                                                reac.subtract_metabolites({'h_c':balance_test['H']})
                                            elif reac_comp == 'p':
                                                reac.subtract_metabolites({'h_p':balance_test['H']})
                                            elif reac_comp == 'e':
                                                reac.subtract_metabolites({'h_e':balance_test['H']})
                                            else:
                                                perform_deletion = False
                                                break
                                        # ..............................
                                        # TODO:
                                        #    fix other possible problems
                                        # ..............................

                                        # finally, check balance again (continue only if fixed, else break)
                                        if len(reac.check_mass_balance()) > 0:
                                            perform_deletion = False
                                            break

                                    else:
                                        continue

                                # if not problems are found, duplicate is removed
                                if perform_deletion:
                                    model = model_del.copy()
                                    print(F'\tDuplicate metabolite {del_meta_id} found. Replaced with {keep_meta.id}.')
                                # if problems are not solvable, duplicate is kept and only reported
                                else:
                                    print(F'\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tAutomated deletion not possible due to problems with consistency.')

                        # else, metabolite is kept
                        #       since it might be an isomer, elongation, or other explanation
                        #       for the same annotation
                        else:
                            print(F'\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tKept, as reaction containing both metabolites was found.')


                # ... or only report duplicates
                else:
                    print(F'\tDuplicate metabolite(s) {", ".join(mnx[1]["id"].tolist())} found.')


    return model


def resolve_duplicates(model:cobra.Model, check_reac:bool=True, 
                       check_meta:Literal['default','exhaustive','skip']='default', 
                       replace_dupl_meta:bool=True, remove_unused_meta:bool=False, 
                       remove_dupl_reac:bool=True) -> cobra.Model:
    """Resolve and remove (optional) duplicate metabolites and reactions in the model.

    Args:
        - model (cobra.Model): 
            The model loaded with COBRApy.
        - check_reac (bool, optional): 
            Whether to check reactions for duplicates. 
            Defaults to True.
        - check_meta (Literal['default','exhaustive','skip'], optional): 
            Whether to check for duplicate metabolites. 
            Defaults to 'default'.
        - replace_dupl_meta (bool, optional): 
            Option to replace/remove duplicate metabolites. 
            Defaults to True.
        - remove_unused_meta (bool, optional): 
            Option to remove unused metabolites. 
            Defaults to False.
        - remove_dupl_reac (bool, optional): 
            Option to combine/remove duplicate reactions. 
            Defaults to True.

    Returns:
        cobra.Model: 
            The (edited) model.
    """

    # resolve duplicate metabolites
    if check_meta == 'default':
        # resolve duplicates starting with the metanetx.chemical database identifiers
        model = resolve_duplicate_metabolites(model, replace=replace_dupl_meta)
    elif check_meta == 'exhaustive':
        # resolve duplicates by starting at every database identifer one after another
        # note: bigg and sbo are skipped as sbo gives not much information and bigg is
        #       usually the one that differs (naming issue)
        anno_types = set()
        # get all database annotation types present in the model
        for m in model.metabolites:
            anno_types = anno_types | set(m.annotation.keys())
        for colname in [_ for _ in anno_types if not _ in ['bigg.metabolite','sbo']]:
            model = resolve_duplicate_metabolites(model,colname,replace=replace_dupl_meta)
    elif check_meta == 'skip':
        print('\tSkip check for duplicate metabolites.')
    else:
         warnings.warn(F'Warning: Unknown option for metabolites duplicate checking {check_meta}. Search for metabolite duplicates skipped.')

    # remove now unused metabolites
    if remove_unused_meta:
        model,removed = cobra.manipulation.delete.prune_unused_metabolites(model)
        print(F'\tThe following metabolites () have been removed: {", ".join([x.id for x in removed])}')

    # resolve duplicate reactions
    if check_reac:
        model = resolve_duplicate_reactions(model, based_on='reaction', remove_reac =remove_dupl_reac)

    return model


