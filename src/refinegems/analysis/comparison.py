#!/usr/bin/env python
"""Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models.
All other stats are shown in the memote report.
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import logging
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from cobra import Model as cobraModel
from libsbml import Model as libModel
from libsbml import BIOLOGICAL_QUALIFIER, BQB_IS, BQB_IS_HOMOLOG_TO
from pathlib import Path
from typing import Literal, Union
from upsetplot import from_memberships, UpSet
from venn import venn

from ..classes.reports import MultiSBOTermReport, SBOTermReport
from ..curation.miriam import get_set_of_curies

################################################################################
# setup logging
################################################################################

logger = logging.getLogger(__name__)

################################################################################
# functions
################################################################################


def sbo_terms(models: list[libModel], rename: Union[list[str],None]=None) -> MultiSBOTermReport:
    """Analyse and compare the SBO term annotations of a given list
    of models.

    Args:
        - models (list[libModel]):
            A list containing models loaded with libSBML.
        - rename (list[str], optional):
            Rename model ids to custom names.
            Defaults to None.

    Returns:
        MultiSBOTermReport:
            A :py:class:`~refinegems.classes.reports.MultiSBOTermReport` instance.
    """

    sboanalyses = []
    for m, name in zip(models, rename):
        sboanalyses.append(SBOTermReport(m, name))

    return MultiSBOTermReport(sboanalyses)

###
# Model entity comparison
###
# @TODO add handling of IDs as they can differ between models for the suffix, e.g. _c vs. __61__c__63__
def get_entity_curie_set_per_db(model: libModel, entity: Literal['genes', 'metabolites', 'reactions', 'pathways'], db: Union[str, None], include_homologs: bool=False) -> Union[set, None]:
    """Get set of CURIEs for one entity type for a specific database

    Args:
        model (libModel): Model laoded with libSBML
        entity (str): String specifying the entity type. One of: metabolites|reactions|pathways
        db (str): Specifies the database the entity identifier should come from. 
        Needs to be one of the keys in bioregistry.get_prefix_map().

    Returns:
        Union[set, None]: Set of identifiers found in the model for the given entity type and database if available, else None.

    Raises:
        - ValueError: If entity type is unknown
    """
    list_of = []

    match entity:
        case 'genes':
           list_of = model.getPlugin('fbc').getListOfGeneProducts()
        case 'metabolites':
            list_of = model.getListOfSpecies()
        case 'reactions':
            list_of = model.getListOfReactions()
        case 'pathways':
            try:
                list_of = model.getPlugin('groups').getListOfGroups()
            except AttributeError:
                logging.warning(f'Model {model.getId()} does not contain any pathway information.')
                return None
        case _:
            raise ValueError(f"Unknown entity: {entity!r}. Expected one of: {', '.join(['genes', 'metabolites', 'reactions', 'pathways'])}")

    entity_curie_set = set()

    for e in list_of:

        if db:
            cvterms = e.getCVTerms()

            for cvt in cvterms:
                # Only get CURIEs from IS relationship or include CURIEs with BQB_IS_HOMOLOG_TO if include_homologs=True
                curie_condition = (
                    (cvt.getQualifierType() == BIOLOGICAL_QUALIFIER) and ((cvt.getBiologicalQualifierType() == BQB_IS) 
                    or (cvt.getBiologicalQualifierType() == BQB_IS_HOMOLOG_TO)) 
                    if include_homologs else 
                    (cvt.getQualifierType() == BIOLOGICAL_QUALIFIER) and (cvt.getBiologicalQualifierType() == BQB_IS)
                    )
                if curie_condition:
                    current_curies = [cvt.getResourceURI(i) for i in range(cvt.getNumResources())]
                    prefix2id = get_set_of_curies(current_curies)[0]
                    for key in prefix2id:
                        if db in key:
                            id_set = prefix2id.get(key)
                            for i in id_set:
                                if i != 'NaN': # Exclude NaN identifiers
                                    entity_curie_set.add(i)
        else:
            entity_curie_set.add(e.getId())

    return entity_curie_set


def plot_venn(
    models: list[cobraModel], entity: str, perc: bool = False, rename=None
) -> matplotlib.axes.Axes:
    """Creates Venn diagram to show the overlap of model entities

    Args:
        - models (list[cobraModel]):
            Models loaded with cobrapy
        - entity (str):
            Compare on metabolite|reaction
        - perc (bool, optional):
            True if percentages should be used.
            Defaults to False.
        - rename (dict, optional):
            Rename model ids to custom names.
            Defaults to None.

    Returns:
        matplotlib.axes.Axes:
            Venn diagram
    """
    intersec = {} # model ID : list of entity IDs
    for model in models:
        reas = [] # List of current entity IDs
        if entity == "metabolite":
            for rea in model.metabolites:
                reas.append(rea.id)
        if entity == "reaction":
            for rea in model.reactions:
                reas.append(rea.id)
        if rename is not None:
            intersec[rename[model.id]] = set(reas)
        else:
            intersec[model.id] = set(reas)
    if perc:
        fig = venn(intersec, fmt="{percentage:.1f}%")
    else:
        fig = venn(intersec)
    return fig


def plot_db_entity_overlap(
    models: list[cobraModel], entity: Literal['genes', 'metabolites', 'reactions', 'pathways'], db: Union[str, None]=None, 
    cmap: Union[list[tuple], None]=None, rename: list[str]=None, include_homologs: bool=False, 
    venn_kwargs: dict={'fmt': "{percentage:.1f}%", 'legend_loc':'lower right'}, 
    upset_min_subset_size: int=15,
    outfile_suffix: str='pdf', outdir: str='./'
) -> matplotlib.axes.Axes:
    """Creates Venn diagram (<= 4 models) or UpSet plot to show the overlap of model entities (based on a specific database)

    Args:
        - models (list[cobraModel]): 
             Models loaded with libSBML
        - entity (Literal[genes, metabolites, reactions, pathways]): 
            Entity to compare on.
            Can be one of [genes, metabolites, reactions, pathways].
        - db (str, optional)
            Specifies the database the entity identifier should come from. 
            Needs to be one of the keys in bioregistry.get_prefix_map().
            Defaults to None, which uses model internal IDs.
        - cmap (list[tuple]): 
            Specify colour map to be used by venn/upsetplot or None to use default colours. 
            Defaults to None.
        - rename (str, optional): 
            Rename model ids to custom names.
            Defaults to None.
        - include_homologs (bool, optional): 
            Specifies if homologs should be included in the comparison. 
            Defaults to False.
        - venn_kwargs (dict, optional): 
            Dictionary containing details for plotting the venn plots. 
            Defaults to {'fmt': "{percentage:.1f}%", 'legend_loc':'lower right'}.
        - upset_min_subset_size (int, optional): 
            Minimum subset size to be shown in the UpSet plot. 
            Defaults to 15.
        - outfile_suffix (str, optional): 
            File extension to be used for saved figures. 
            Defaults to 'pdf'.
        - outdir (str, optional): 
            Directory to save resulting figure to. 
            Defaults to './'.

    Returns:
        matplotlib.axes.Axes: 
            Figure resulting from venn/upsetplot
    """

    # Set-up filepath and folder to store resulting figure in
    if outdir:
        if isinstance(outdir, str): outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        filename = f'{entity}_{db}_comparison_plot' if db else f'{entity}_ID_comparison_plot'
        filepath = Path(outdir, f'{filename}.{outfile_suffix}')

    # If no renaming list is given, use model IDs
    if not rename: rename = [model.id for model in models]
    
    # Get CURIE set for specifed entity and database for each model
    model2ids = {model_id2use: get_entity_curie_set_per_db(model, entity, db, include_homologs)
                 for model, model_id2use in zip(models, rename)} # model ID : set of entity IDs

    # Filter out Nones
    model2remove = []
    for m in model2ids.keys():
        if model2ids[m] is None:
            comp_method = f'{db} {entity}' if db else f'{entity}'
            logging.warning(f'Skipping model {m} for {comp_method} comparison. No {entity} information available in the model.')
            model2remove.append(m)
    model2ids = {m: model2ids[m] for m in model2ids.keys() if m not in model2remove}

    # Generate and save plot
    if len(models) <= 4:
        fig = venn(model2ids, cmap=cmap, **venn_kwargs)
    else:
        # Prepare data for UpSet plot
        all_ids = sorted(set.union(*model2ids.values()))
        memberships = []
        for id_ in all_ids:
            present_in = [mid for mid, ids in model2ids.items() if id_ in ids]
            memberships.append(present_in)

        upset_data = from_memberships(memberships, data=[1]*len(memberships))
        fig = plt.figure()
        upset = UpSet(upset_data, subset_size='count', min_subset_size=upset_min_subset_size, show_counts=True )

        # Set colours per degree if provided
        if cmap:
            # Determine the maximum degree present in the intersections
            max_degree = max(len(idx) for idx in upset_data.index)
            
            # Colour bars per degree
            for degree in range(1, max_degree + 1):
                colour = cmap[(degree - 1) % len(cmap)] # Cycle through provided colours if not enough
                upset.style_subsets(min_degree=degree, facecolor=colour)

        # Plot UpSetPlot
        upset.plot(fig=fig)
    if outdir:
        plt.savefig(filepath, dpi=300)
        plt.close()
    else:
        logging.info(f'No ouput directory given. Resulting figure returned but NOT saved.')

    # Also return figure
    return fig
