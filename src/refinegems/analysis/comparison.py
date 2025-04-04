#!/usr/bin/env python
"""Provides functions to compare and visualize multiple models

Can mainly be used to compare growth behaviour of multiple models.
All other stats are shown in the memote report.
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import matplotlib
import matplotlib.colors

from cobra import Model as cobraModel
from libsbml import Model as libModel
from venn import venn

from ..classes.reports import MultiSBOTermReport, SBOTermReport

################################################################################
# functions
################################################################################


def sbo_terms(models: list[libModel]) -> MultiSBOTermReport:
    """Analyse and compare the SBO term annotations of a given list
    of models.

    Args:
        - models (list[libModel]):
            A list containing models loaded with libSBML.

    Returns:
        MultiSBOTermReport:
            A :py:class:`~refinegems.classes.reports.MultiSBOTermReport` instance.
    """

    sboanalyses = []
    for m in models:
        sboanalyses.append(SBOTermReport(m))

    return MultiSBOTermReport(sboanalyses)


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
    intersec = {}
    for model in models:
        reas = []
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
