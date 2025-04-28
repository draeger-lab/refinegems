"""Classes to generate, handle, manipulate and save reports."""

__author__ = "Carolin Brune, Famke Baeuerle, Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################

import cobra
import copy
import logging
import math
import matplotlib
import matplotlib.figure
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
import numpy as np
import pandas as pd
import re
import seaborn as sns
import warnings

from importlib.resources import files
from itertools import chain
from libsbml import Model as libModel
from pathlib import Path
from typing import Literal, Union

from ..analysis.investigate import (
    get_mass_charge_unbalanced,
    get_orphans_deadends_disconnected,
    get_reac_with_gpr,
    get_reactions_per_sbo,
)
from ..utility.util import test_biomass_presence
from ..developement.decorators import *
from ..utility.io import search_sbo_label

################################################################################
# variables
################################################################################

KEGG_GLOBAL_PATHWAY = {
    "01100": "Metabolic pathways",
    "01110": "Biosynthesis of secondary metabolites",
    "01120": "Microbial metabolism in diverse environments",
}  #: :meta:

KEGG_OVERVIEW_PATHWAY = {
    "01200": "Carbon metabolism",
    "01210": "2-Oxocarboxylic acid metabolism",
    "01212": "Fatty acid metabolism",
    "01230": "Biosynthesis of amino acids",
    "01232": "Nucleotide metabolism",
    "01250": "Biosynthesis of nucleotide sugars",
    "01240": "Biosynthesis of cofactors",
    "01220": "Degradation of aromatic compounds",
}  #: :meta:

KEGG_METABOLISM_PATHWAY = files("refinegems.data.pathway").joinpath(
    "KEGG_pathway_metabolism.csv"
)  #: :meta hide-value:
KEGG_METABOLISM_PATHWAY_DATE = "6. July 2023"  #: :meta:

################################################################################
# classes
################################################################################


class Report:
    pass


class SingleGrowthSimulationReport(Report):
    """Report for a single growth simulation, one media against one model.

    Attributes:
        - model_name:
            Name of the model.
        - medium_name:
            Name of the medium.
        - growth_value:
            Simulated growth value.
        - doubling_time:
            Simulated doubling time.
        - additives:
            List of substances, that were added.
        - no_exchange:
            List of substances that normally would be found in the media
            but have been removed, as they are not part of the exchange reactions
            of the model.
    """

    def __init__(
        self,
        model_name=None,
        medium_name=None,
        growth_value=None,
        doubling_time=None,
        additives=None,
        no_exchange=None,
    ):
        self.model_name = model_name
        self.medium_name = medium_name
        self.growth_value = growth_value
        self.doubling_time = doubling_time
        self.additives = additives
        self.no_exchange = no_exchange

    def __str__(self):
        return (
            f"model: {self.model_name}\n"
            f"medium: {self.medium_name}\n"
            f"growth: {self.growth_value}\n"
            f"doubling time: {self.doubling_time}\n"
            f"additives: {self.additives}\n"
            f"no_exchange: {self.no_exchange}\n"
        )

    def to_dict(self) -> dict:
        """Transform the information into a dictionary.

        Returns:
            dict:
                The information of the report as a dictionary.
        """

        return {
            "model_name": self.model_name,
            "medium_name": self.medium_name,
            "growth_value": self.growth_value,
            "doubling_time": self.doubling_time,
            "additives": self.additives if len(self.additives) > 0 else None,
            "no_exchange": self.no_exchange,
        }


class GrowthSimulationReport(Report):
    """Report for the growth simulation analysis.

    Attributes:
        - reports:
            List of the report for the single growth analysis.
        - model:
            List of the model names.
        - media:
            List of the media names.
    """

    def __init__(self, reports: list[SingleGrowthSimulationReport] = None):

        self.reports = reports if reports else []
        self.models = set([_.model_name for _ in reports]) if reports else set()
        self.media = set([_.medium_name for _ in reports]) if reports else set()

    def __str__(self) -> str:

        return "\n\n".join(str(_) for _ in self.reports)

    def add_sim_results(self, new_rep: SingleGrowthSimulationReport):
        """Add a new single growth report to the reports list

        Args:
            - new_rep (SingleGrowthSimulationReport):
                The new simulation report.
        """

        self.reports.append(new_rep)
        self.models.add(new_rep.model_name)
        self.media.add(new_rep.medium_name)

    def to_table(self) -> pd.DataFrame:
        """Return a table of the contents of the report.

        Returns:
            pd.DataFrame:
                The table containing the information in the report.
        """

        l = []
        for report in self.reports:
            l.append(report.to_dict())
        return pd.DataFrame(l)

    def plot_growth(
        self, unit: Literal["h", "dt"] = "dt", color_palette: str = "YlGn"
    ) -> matplotlib.figure.Figure:
        """Visualise the contents of the report.

        .. note::

            Please keep in mind that the figure does not show unrealistically high and minicules values to zero.
            However, all values are contained within the table one can get via
            :py:func:`~refinegems.classes.reports.GrowthSimulationReport.to_table`.

        Args:
            - unit (Literal['h','dt'], optional):
                Set the unit to plot.
                Can be doubling time in minutes ('dt') or growth rates in mmol/gDWh ('h').
                Defaults to 'dt'.
            - color_palette (str, optional):
                A colour gradient from the matplotlib library.
                If the name does not exist, uses the default.
                Defaults to 'YlGn'.

        Returns:
            matplotlib.figure.Figure:
                The plotted figure.
        """

        def plot_growth_bar(
            xdata: list[str],
            xlab: str,
            ydata: list[float],
            ylab: str,
            title: str,
            color_palette: str = "YlGn",
        ) -> matplotlib.figure.Figure:
            """Helper function to plot the bar plot for the growth visualisation.

            Args:
                - xdata (list[str]):
                    List of the x-axis data (medium or model names).
                - xlab (str):
                    The x-axis label.
                - ydata (list[float]):
                    List of thr y-axis data (the values).
                - ylab (str):
                    The y-axis label.
                - title (str):
                    The title of the plot.
                - color_palette (str, optional):
                    A colour gradient from the matplotlib library.
                    If the name does not exist, uses the default.
                    Defaults to 'YlGn'.

            Returns:
                matplotlib.figure.Figure:
                    The plotted figure.
            """

            # create colour gradient
            try:
                cmap = matplotlib.colormaps[color_palette]
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.colormaps["YlGn"]

            # set up the figure
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])

            # clean-up data
            ydata = [_ if _ > 0.0 else 0.0 for _ in ydata]
            ydata = [_ if _ < 1000.0 else 0.0 for _ in ydata]

            # construct the plot
            max_ydata = max(ydata)
            if max_ydata <= 0:
                warnings.warn(
                    "Model is not able to grow sensible on any medium. Returning empty figure."
                )
                return fig
            cont = ax.bar(
                list(xdata), ydata, color=cmap([_ / max_ydata for _ in ydata])
            )

            # set labels and others
            ax.bar_label(cont, fmt="%.2f", color="black", padding=1.0)
            ax.set_ylabel(ylab, labelpad=12)
            ax.set_xlabel(xlab, labelpad=12)
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0], xlims[1] + (xlims[1] * 0.03))
            ax.set_title(title)

            return fig

        def plot_growth_heatmap(
            data: pd.DataFrame, color_palette: str = "YlGn"
        ) -> matplotlib.figure.Figure:
            """Helper function to plot the heatmap for the growth visualisation.

            Args:
                - data (pd.DataFrame):
                    The table containing the data to be plotted.
                    Needs to have the columns 'medium', 'model' and one for the growth values.
                - color_palette (str, optional):
                    A colour gradient from the matplotlib library.
                    If the name does not exist, uses the default.
                    Defaults to 'YlGn'.

            Returns:
                matplotlib.figure.Figure:
                    The plotted figure.
            """

            # clean up + transform data
            growth = data.set_index(["medium", "model"]).sort_index().T.stack()
            growth.columns.name = None
            growth.index.names = (None, None)
            growth.index.name = None
            growth.index = growth.index.get_level_values(1)

            # over / under (meaningful) values
            growth[growth > 1000] = 0
            growth[growth < 0] = 0
            growth.replace([np.inf, -np.inf], 0, inplace=True)
            over_growth = growth.max().max() + 6
            growth.replace(np.nan, over_growth, inplace=True)
            under_growth = growth.min().min() - 5
            vmin = (
                under_growth if under_growth > 1e-5 else 1e-5
            )  # Use same threshhold as in find_missing_essential in growth
            vmax = over_growth - 1

            # annotations
            annot = growth.copy()
            annot = annot.round().astype(int)
            annot[annot < 1e-5] = ""
            annot.replace(over_growth.round().astype(int), "", inplace=True)

            # setting the colours
            try:
                cmap = matplotlib.colormaps[color_palette]
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.colormaps["YlGn"]
            cmap.set_under("white")  # too low / no growth
            cmap.set_over("white")  # no data

            # plot the heatmap
            fig, ax = plt.subplots(figsize=(10, 8))

            zm = np.ma.masked_where(growth.T != over_growth, growth.T)
            x = np.arange(len(growth.T.columns) + 1)
            y = np.arange(len(growth.T.index) + 1)

            res = sns.heatmap(
                growth.T,
                annot=annot.T,
                annot_kws={"fontsize": 15},
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
                linewidth=0.5,
                cbar_kws={
                    "orientation": "vertical",
                    "label": "Doubling time [min]",
                    "extend": "min",
                    "extendrect": True,
                },
                ax=ax,
                fmt="",
            )

            # labels
            rotation = 40 if len(growth.index) > 3 else 0
            plt.tick_params(
                rotation=0, bottom=False, top=False, left=False, right=False
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha="right")
            ax.set_ylabel("Media", fontsize=15)
            ax.set_xlabel("Models", fontsize=15)

            # hatches for no-data-values
            plt.pcolor(x, y, zm, hatch="x", alpha=0.0)

            # drawing a frame
            for _, spine in res.spines.items():
                spine.set_visible(True)
                spine.set_linewidth(1)
                spine.set_edgecolor("grey")

            # extra legend
            handles = []
            handles.append(
                mpatches.Rectangle(
                    (0, 0), 0, 0, color="white", ec="grey", label="No growth"
                )
            )
            handles.append(
                mpatches.Rectangle(
                    (0, 0), 0, 0, color="white", ec="grey", hatch="xxx", label="No data"
                )
            )
            fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.9, 0.05))

            return fig

        # match the unit
        match unit:
            case "h":
                unit_text = r"growth rate $[\frac{mmol}{gDWh}]$"
            case "dt":
                unit_text = "doubling time [min]"
            case _:
                raise ValueError(f"Unknown input for unit: {unit}")

        # multiple models vs one medium
        if len(self.models) > 1 and len(self.media) == 1:

            # collect data
            xdata = [_.model_name for _ in self.reports]
            xlab = "models"
            ydata = [
                _.growth_value if unit == "h" else _.doubling_time for _ in self.reports
            ]
            ylab = unit_text
            title = (
                f"Growth simulation on {next(iter(self.media))} for different models"
            )

            # plot
            return plot_growth_bar(xdata, xlab, ydata, ylab, title, color_palette)

        # one model vs mutiple media
        elif len(self.models) == 1 and len(self.media) > 1:

            # collect data
            xdata = [_.medium_name for _ in self.reports]
            xlab = "media"
            ydata = [
                _.growth_value if unit == "h" else _.doubling_time for _ in self.reports
            ]
            ylab = unit_text
            title = (
                f"Growth simulation for {next(iter(self.models))} on different media"
            )

            # plot
            return plot_growth_bar(xdata, xlab, ydata, ylab, title, color_palette)

        # one medium, one model case - just to make it usable for all inputs
        elif len(self.models) == 1 and len(self.media) == 1:

            # collect data
            xdata = [_.medium_name for _ in self.reports]
            xlab = "medium"
            ydata = [
                _.growth_value if unit == "h" else _.doubling_time for _ in self.reports
            ]
            ylab = unit_text
            title = (
                f"Growth simulation for {next(iter(self.models))} on medium {xdata[0]}"
            )

            # plot
            return plot_growth_bar(xdata, xlab, ydata, ylab, title, color_palette)

        # multiple vs multiple
        elif len(self.models) > 1 and len(self.media) > 1:

            data = pd.DataFrame(
                {
                    "model": [_.model_name for _ in self.reports],
                    "medium": [_.medium_name for _ in self.reports],
                    "growth": [
                        _.growth_value if unit == "h" else _.doubling_time
                        for _ in self.reports
                    ],
                }
            )

            return plot_growth_heatmap(data, color_palette)

        # problematic case
        else:
            raise IndexError(
                "Can only plot growth if at least one model and one medium are present."
            )

    def save(
        self,
        to: str,
        how: Literal["dir"] = "dir",
        check_overwrite: bool = True,
        color_palette: str = "YlGn",
    ):
        """Save the report.

        Current options include:

        - 'dir': save the report to a directory, including a txt and two graphics

        Args:
            - to (str):
                Path to a directory to save the report to.
            - how (Literal['dir'], optional):
                How to save the report.
                For options see functions description.
                Defaults to 'dir'.
            - check_overwrite (bool, optional):
                Flag to choose to check for existing directory/files of same name
                or just to overwrite them. Defaults to True.
            - color_palette (str, optional):
                A colour gradient from the matplotlib library.
                If the name does not exist, uses the default.
                Defaults to 'YlGn'.

        Raises:
            - ValueError: If the parameter 'how' is given something unexpected.
        """

        match how:
            # save to a new directory
            case "dir":
                # create directory to save report to
                dir_path = Path(to, "GrowthSimReport")
                dir_path.mkdir(parents=True, exist_ok=check_overwrite)
                # save the report
                with open(Path(dir_path, "report.txt"), "w") as f:
                    f.write(str(self))
                # save visualisation for doubling time
                fig_dt = self.plot_growth(color_palette=color_palette)
                fig_dt.savefig(Path(dir_path, "report_vis_dt.png"), bbox_inches="tight")
                # save visualisation for growth rate
                fig_dt = self.plot_growth(unit="h", color_palette=color_palette)
                fig_dt.savefig(Path(dir_path, "report_vis_h.png"), bbox_inches="tight")

            case _:
                raise ValueError(
                    f'Unknow input for parameter "how": {how}.\n Cannot save report. Abort.'
                )


class KEGGPathwayAnalysisReport(Report):
    """Report for the KEGG pathway analysis.

    Attributes:
        - total_reac:
            An integer for the total number of reactions in the model.
        - kegg_count:
            An integer as a counter for the KEGG pathway annotations.
        - kegg_global:
            Dictionary of global KEGG IDs and their counts.
        - kegg_over:
            Dictionary of overvire KEGG IDs and their counts.
        - kegg_rest:
            Dictionary of the remaining KEGG IDs and their counts.
    """

    def __init__(
        self,
        total_reac=None,
        kegg_count=None,
        kegg_global=None,
        kegg_over=None,
        kegg_rest=None,
    ) -> None:

        # super().__init__()
        # general counts
        self.total_reac = total_reac
        self.kegg_count = kegg_count

        # kegg pathways
        self.kegg_global = kegg_global
        self.kegg_over = kegg_over
        self.kegg_paths = kegg_rest

    def visualise_kegg_counts(
        self, colors: list[str] = ["lightgreen", "darkgreen"]
    ) -> plt.figure:
        """Visualise the amounts of reaction with and without
        KEGG pathway annotation.

        Args:
            - colors (list[str], optional):
                List of two colours used for the plotting.
                If wrong number or non-matplotlib colours are given, sets its to the default.
                Defaults to 'lightgreen' and 'darkgreen'.

        Returns:
            plt.figure:
                The resulting plot.
        """

        # validate colors
        if (
            len(colors) != 2
            or not matplotlib.colors.is_color_like(colors[0])
            or not matplotlib.colors.is_color_like(colors[1])
        ):
            warnings.warn(
                "Unknown colors or false amount of colors for pie chart. Resume using default values."
            )
            colors = ["lightgreen", "lightskyblue"]

        # generate the plot
        explode = (0.0, 0.1)
        fig, ax = plt.subplots()
        values = [self.kegg_count, self.total_reac - self.kegg_count]
        labels = ["yes", "no"]
        ax.pie(
            values,
            autopct=lambda pct: "{:.1f}%\n({:.0f})".format(
                pct, (pct / 100) * sum(values)
            ),
            colors=colors,
            explode=explode,
            shadow=True,
            startangle=90,
        )
        ax.legend(labels, title="KEGG\npathway")

        return fig

    def visualise_kegg_pathway(
        self,
        plot_type: Literal["global", "overview", "high", "existing"] = "global",
        label: Literal["id", "name"] = "id",
        color_palette: str = "YlGn",
    ) -> plt.figure:
        """Visualise the KEGG pathway identifiers present.

        Depending on the :plot_type:, different levels of pathway identifiers
        are plotted:

        - global: check and plot only the global pathway identifiers
        - overview: check and plot only the overview pathway identifiers
        - high: check and plot all identifiers grouped by their high level pathway identifiers. This option uses label=name, independedly of the input
        - all: check and plot all identifiers

        Args:
            - plot_type (Literal["global","overview","high","existing"], optional):
                Type of plot, explaination see above. Defaults to 'global'.
            - label (Literal["id","name"], optional):
                Type of the label. If 'id', uses the KEGG pathway IDs,
                if 'name', uses the pathway names. Defaults to 'id'.
            - color_palette (str, optional):
                A colour gradient from the matplotlib library.
                If the name does not exist, uses the default.
                Defaults to 'YlGn'.

        Returns:
            plt.figure:
                The plotted visualisation.
        """

        # get data and KEGG pathway label mapping
        # for the given plot type
        match plot_type:
            case "global":
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = "global identifiers"
            case "overview":
                data = self.kegg_over
                label_map = KEGG_OVERVIEW_PATHWAY
                title_type = "overview identifiers"
            case "high":
                label = "name"
                data = self.kegg_paths
                label_map = (
                    pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str)
                    .set_index("id")[["group"]]
                    .to_dict()["group"]
                )
                title_type = "grouped identifiers"
            case "existing":
                data = self.kegg_paths
                label_map = (
                    pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str)
                    .set_index("id")[["specific"]]
                    .to_dict()["specific"]
                )
                title_type = "identifiers in model"
            case _:
                warnings.warn(
                    f'Unknown option for plot_type, choosing "global" istead: {plot_type}'
                )
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = "global identifiers"

        # create the label
        match label:
            case "id":
                for k in label_map:
                    if k not in data and plot_type != "existing":
                        data[k] = 0
            case "name":
                old_data = data
                data = {}
                for k in label_map:
                    if k not in old_data:
                        if plot_type != "existing" and label_map[k] not in data:
                            data[label_map[k]] = 0
                    else:
                        if label_map[k] not in data:
                            data[label_map[k]] = old_data[k]
                        else:
                            data[label_map[k]] += old_data[k]

            case _:
                warnings.warn(f'Unknown input for label: {label}. Using "id" instead.')
                for k in label_map:
                    if k not in data:
                        data[k] = 0

        data = pd.DataFrame(data.items(), columns=["label", "counts"]).sort_values(
            "counts"
        )
        cdata = data.counts.values
        ldata = data.label.values

        # create the graph
        # ----------------
        # setting the colours
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]
        # set up the figure
        fig = plt.figure()
        if "existing" == plot_type:
            ax = fig.add_axes([0, 0, 5, 6])
            # construct the plot
            cont = ax.barh(
                ldata, cdata, color=cmap([x / max(cdata) for x in cdata]), label=ldata
            )
            ax.bar_label(cont, fmt="%d", color="black", padding=1.0)
            ax.set_ylabel("KEGG pathway")
            ax.set_xlabel("Number of reaction annotations")
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0], xlims[1] + (xlims[1] * 0.03))
            ax.set_title(f"Pathway analysis with KEGG: {title_type}")
        else:

            ax = fig.add_axes([0, 0, 1, 1])
            # construct the plot
            cont = ax.barh(
                ldata, cdata, color=cmap([x / max(cdata) for x in cdata]), label=ldata
            )
            ax.bar_label(cont, fmt="%d", color="black", padding=1.0)
            ax.set_ylabel("KEGG pathway", labelpad=12)
            ax.set_xlabel("Number of reaction annotations", labelpad=12)
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0], xlims[1] + (xlims[1] * 0.03))
            ax.set_title(f"Pathway analysis with KEGG: {title_type}")

        return fig

    def save(self, dir: str, colors: str = "YlGn") -> None:
        """Save the content of the report as plots.

        Args:
            - dir (str):
                Path to a directory to save the output directory with all the plot in.
            - colors(str,optional):
                Colour palette for the plots.
                Should be a valid name of a matplotlib sequential colour palette.
        """

        # collect all produced file in one directory
        try:
            Path(dir, "pathway-analysis").mkdir(parents=True, exist_ok=False)
            print(f'Creating new directory {str(Path(dir,"pathway-analysis"))}')
        except FileExistsError:
            print("Given directory already has required structure.")

        # create and save plots
        # a) for the counts
        if self.total_reac and self.kegg_count:
            count_fig = self.visualise_kegg_counts()
            count_fig.savefig(
                str(Path(dir, "pathway-analysis", "kegg_anno_counts.png")),
                bbox_inches="tight",
            )
        # b) for the actual pathways
        # 1.) global KEGG IDs
        if self.kegg_global:
            # with id
            fig = self.visualise_kegg_pathway(
                plot_type="global", label="id", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_global_id.png")),
                bbox_inches="tight",
            )
            # with name
            fig = self.visualise_kegg_pathway(
                plot_type="global", label="name", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_global_name.png")),
                bbox_inches="tight",
            )
        # 2.) Overview KEGG IDs
        if self.kegg_over:
            # with id
            fig = self.visualise_kegg_pathway(
                plot_type="overview", label="id", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_overview_id.png")),
                bbox_inches="tight",
            )
            # with name
            fig = self.visualise_kegg_pathway(
                plot_type="overview", label="name", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_overview_name.png")),
                bbox_inches="tight",
            )
        # 3.) rest
        if self.kegg_paths:
            # grouped by high-level terms
            fig = self.visualise_kegg_pathway(
                plot_type="high", label="name", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_high.png")),
                bbox_inches="tight",
            )
            # all with id
            fig = self.visualise_kegg_pathway(
                plot_type="existing", label="id", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_existing_id.png")),
                bbox_inches="tight",
            )
            # all with name
            fig = self.visualise_kegg_pathway(
                plot_type="existing", label="name", color_palette=colors
            )
            fig.savefig(
                str(Path(dir, "pathway-analysis", "pathway_existing_name.png")),
                bbox_inches="tight",
            )


class AuxotrophySimulationReport(Report):
    """Report for the auxotrophy simulation.

    Attributes:
        - simulation_results:
            The data of the simulation.
    """

    def __init__(self, results) -> None:
        # super().__init__()
        self.simulation_results = results

    # auxotrophy sim visualisation
    def visualise_auxotrophies(
        self, color_palette: str = "YlGn", save: Union[None, str] = None
    ) -> Union[None, matplotlib.figure.Figure]:
        """Visualise and/or save the results of the :py:func:`~refinegems.analysis.growth.test_auxotrophies` function.

        Args:
            - res (pd.DataFrame):
                The output of  :py:func:`~refinegems.analysis.growth.test_auxotrophies`.
            - color_palette (str, optional):
                A name of a seaborn gradient color palette.
                In case name is unknown, takes the default. Defaults to 'YlGn'.
            - save (None | str, optional):
                Path to a directory, if the output shall be saved. Defaults to None (returns the figure).

        Returns:
            (1) Case: ``save = str``
                    None: No
                        return, as the visulaisation is directly saved.

            (2) Case: ``save = None``
                    matplotlib.figure.Figure:
                        The plotted figure.
        """

        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        # set up the figure
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])

        # create heatmap
        sns.heatmap(
            self.simulation_results,
            ax=ax,
            cmap=cmap,
            cbar_kws={"label": "flux"},
            annot=True,
            fmt=".2f",
        )

        # add labels
        ax.set_ylabel("amino acid", labelpad=12)
        ax.set_xlabel("medium", labelpad=12)
        ax.set_title("Fluxes for auxotrophy tests")

        # save or return
        if save:
            # save the visualisation of the growth rates
            fig.savefig(Path(save, "auxotrophies_vis.png"), bbox_inches="tight")
        else:
            return fig

    def save(self, dir: str, color_palette: str = "YlGn"):
        """Save the report to a given dictionary.

        Args:
            - dir (str):
                Path to a dictionary.
            - color_palette (str, optional):
                Name of a matplotlib colour palette. Defaults to 'YnGr'.
        """

        # save the visualisation of the growth rates
        self.visualise_auxotrophies(color_palette, save=dir)

        # save the growth rates as tabular information
        self.simulation_results.to_csv(
            Path(dir, "auxotrophies_table.tsv"), sep="\t", index=True
        )


class SourceTestReport(Report):
    """Report for the source test (:py:func:rg.growth.test_growth_with_source).

    Attributes:
        - results:
            A pd.DataFrame with the results (substances and growth values).
        - element:
            The element the test was performed for.
        - model_name:
            The name of the model, that was tested.
    """

    def __init__(
        self, results: pd.DataFrame = None, element: str = None, model_name: str = None
    ):
        # super().__init__()
        self.results = results
        self.element = element
        self.model_name = model_name

    def visualise(
        self, width: int = 12, color_palette: str = "YlGn"
    ) -> tuple[matplotlib.figure.Figure, pd.DataFrame]:
        """Visualise the results of the source test as a heatmap

        Args:
            - width (int, optional):
                Number of columns to display for the heatmap.
                Number of row is calculated accordingly to fit all values.
                Defaults to 12.
            - color_palette (str, optional):
                Color palette (gradient) for the plot.
                Defaults to 'YlGn'.

        Returns:
            tuple(matplotlib.Figure, pd.DataFrame):
                The heatmap and the legend explaining the heatmap.
        """

        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        cmap.set_under("black")  # too low / no growth
        cmap.set_over("white")  # no data

        # get size of heatmap
        height = math.ceil(len(self.results) / width)
        total_cells = width * height

        # create table for plotting
        data_to_plot = copy.deepcopy(self.results)
        if len(self.results) < total_cells:
            temp = pd.DataFrame.from_records(
                [["empty", None]] * (total_cells - len(self.results)),
                columns=["substance", "growth value"],
            )
            data_to_plot = pd.concat([data_to_plot, temp], ignore_index=True)
        data_to_plot["row"] = list(range(1, height + 1)) * width
        data_to_plot["column"] = list(
            chain.from_iterable([[x] * height for x in range(1, width + 1)])
        )

        # remove unplottable entries
        data_to_plot["growth value"].replace([np.inf, -np.inf], 0, inplace=True)
        over_growth = (
            data_to_plot["growth value"].max()
            + 0.1 * data_to_plot["growth value"].max()
        )
        data_to_plot["growth value"].replace(np.nan, over_growth, inplace=True)
        vmin = 1e-5  # Use same threshhold as in find_missing_essential in growth
        vmax = over_growth - 0.05 * data_to_plot["growth value"].max()

        # set annotations
        annot = data_to_plot.copy()
        annot["growth value"] = annot["growth value"].round(2)
        annot["growth value"] = annot["growth value"].apply(
            lambda x: "" if x < 1e-5 else (x if x < vmax else "X")
        )

        detected_growth = len(
            annot[(annot["growth value"] != "") & (annot["growth value"] != "X")]
        )

        annot = annot.pivot(index="row", columns="column", values="growth value")
        legend = data_to_plot.pivot(index="row", columns="column", values="substance")

        # plot
        ax = sns.heatmap(
            data_to_plot.pivot(index="row", columns="column", values="growth value"),
            linewidth=0.5,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            annot=annot,
            fmt="",
            cbar_kws={"label": r"growth rate $[\frac{mmol}{gDWh}]$"},
        )

        ax.set(xlabel="column", ylabel="row")
        ax.set_title(f"Growth detected with {detected_growth} sources", fontsize=10)
        plt.suptitle(
            f"{self.element}-source growth simulation on model {self.model_name}"
        )

        return (ax.get_figure(), legend)

    def save(self, dir: str, width: int = 12, color_palette: str = "YlGn") -> None:
        """Save the results of the source test.

        Args:
            - dir (str):
                Path to a directory to save the results to.
            - width (int, optional):
                Number of columns for the heatmap.
                Defaults to 12.
            - color_palette (str, optional):
                Color palette (gradient) for the plot.
                Defaults to 'YlGn'.
        """

        # save the list
        self.results.to_csv(
            Path(dir, "source_test_results.csv"), sep=";", header=True, index=False
        )

        # save the visualisation
        fig, leg = self.visualise(width=width, color_palette=color_palette)
        fig.savefig(Path(dir, "source_test_hm.png"), bbox_inches="tight", dpi=400)
        leg.to_csv(
            Path(dir, "source_test_hm_legend.csv"), sep=";", header=True, index=True
        )


class CorePanAnalysisReport(Report):
    """Report for the core-pan analysis.

    Summarises the information and provides functions
    for visualisation.

    Attributes:
        - model:
            The model the report is based on.
        - core_reac:
            List of reactions considered "core".
        - pan_reac:
            List of reactions considered "pan".
        - novel_reac:
            List of reactions considered "novel".
    """

    def __init__(
        self,
        model: cobra.Model,
        core_reac: list[str] = None,
        pan_reac: list[str] = None,
        novel_reac: list[str] = None,
    ):

        # super().__init__()
        # general attributes
        self.model = model
        # reaction attributes
        self.core_reac = core_reac
        self.pan_reac = pan_reac
        self.novel_reac = novel_reac
        # ...

    def get_reac_counts(self):
        """Return a dictionary of the counts of the reactions types (core, pan, novel)."""

        counts = {}
        if self.core_reac:
            counts["core"] = len(self.core_reac)
        else:
            counts["core"] = 0
        if self.pan_reac:
            counts["pan"] = len(self.pan_reac)
        else:
            counts["pan"] = 0
        if self.novel_reac:
            counts["novel"] = len(self.novel_reac)
        else:
            counts["novel"] = 0

        return counts

    def isValid(self, check="reaction-count") -> bool:
        """Check if a certain part of the analysis is valid.

        Currently possible checks:

        - reaction-count :
            check if the number of reactions in the model
            equal the sum of the novel, pan and core reactions

        Args:
            - check (str, optional):
                Describes which part to check. Options are listed above.
                Defaults to 'reaction-count'.

        Raises:
            - ValueError: Unknown string for parameter check.

        Returns:
            bool:
                Result of the check.
        """

        match check:
            case "reaction-count":
                pc_total = sum(self.get_reac_counts().values())
                diff = len(self.model.reactions) - pc_total
                if diff != 0:
                    return False
                else:
                    return True
            case _:
                raise ValueError("Unknown string for parameter check: ", check)

    def visualise_reactions(self) -> matplotlib.figure:
        """Visualise the results of the pan-core analysis for the reactions as a donut chart.

        Returns:
            matplotlib.figure:
                The plot.
        """

        # check the counts
        # ----------------
        counts = self.get_reac_counts()

        if self.isValid("reaction-count"):
            vis_data = list(counts.values())
            vis_label = list(counts.keys())
            vis_color = ["lightgreen", "lightskyblue", "lightcoral"]

        else:
            warnings.warn(
                "Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected."
            )
            vis_data = list(counts.values()).append(
                len(self.model.reactions) - sum(counts.values)
            )
            vis_label = list(counts.keys()).append("discrepancies")
            vis_color = ["lightgreen", "lightskyblue", "lightcoral", "gold"]

        # plot a donut chart
        # ------------------
        fig, ax = plt.subplots()
        wedges, texts, autotexts = ax.pie(
            vis_data,
            autopct=lambda pct: "{:.1f}%\n({:.0f})".format(
                pct, (pct / 100) * sum(vis_data)
            ),
            pctdistance=0.8,
            labeldistance=1.25,
            colors=vis_color,
            radius=1,
            wedgeprops=dict(width=0.4, edgecolor="w"),
        )

        ax.legend(
            wedges,
            vis_label,
            title="Classification",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
        )

        ax.set_title(
            f"Results of core-pan analysis\nof the reactions of model {self.model.id}"
        )

        return fig

    def save(self, dir: str):
        """Save the results inside a PanCoreAnalysisReport object.

        The function creates a new folder 'pan-core-analysis'
        inside the given directory and creates the following documents:

        - table_reactions.tsv : reactions ID mapped to their labels
        - visualise_reactions : donut chart of the values above

        Args:
            - dir (str):
                Path to a directory to save the output to.
        """

        # collect all produced file in one directory
        try:
            Path(dir, "pan-core-analysis/").mkdir(parents=True, exist_ok=False)
            print(f'Creating new directory {str(Path(dir,"pan-core-analysis/"))}')
        except FileExistsError:
            print("Given directory already has required structure.")

        # save the reactions visualisation
        reac_vis = self.visualise_reactions()
        reac_vis.savefig(
            Path(dir, "pan-core-analysis/visualise_reactions.png"), dpi=reac_vis.dpi
        )

        # save table of reactions mapped to characterisation
        if not self.isValid(check="reaction-count"):
            warnings.warn(
                "Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected. Only labbeld reactions will be written to table."
            )
        reac_tab = pd.DataFrame(
            {
                "reaction_id": self.core_reac + self.pan_reac + self.novel_reac,
                "pan-core": (["core"] * len(self.core_reac))
                + (["pan"] * len(self.pan_reac))
                + (["core"] * len(self.novel_reac)),
            }
        )
        reac_tab.to_csv(
            Path(dir, "pan-core-analysis/table_reactions.tsv"), sep="\t", index=False
        )


class ModelInfoReport(Report):
    """Report about the basic information of a given model.

    Attributes:
        - name:
            A string for the name of the model.
        - reac:
            List of the reactions in the model.
        - meta:
            List of the metabolites in the model.
        - gene:
            An int that describes the number of genes in the model.
        - orphans:
            List of metabolite IDs that are considered orphans.
        - deadends:
            List of metabolite IDs that are considered dead-ends.
        - disconnects:
            List of metabolite IDs that are disconnected in the model.
        - mass_charge_unbalanced:
            List of reaction IDs that are unbalanced regarding their mass and their charges.
        - mass_unbalanced:
            List of reaction IDs that are unbalanced regarding their mass only.
        - charge_unbalanced:
            List of reactions IDs that are unbalanced regarding their charges only.
        - pseudo:
            List of pseudoreaction IDs (sinks, demands, exchanges) in the model.
        - normal_with_gpr:
            List of reactions IDs that are normal reactions with gpr.
        - pseudo_with_gpr:
            List of reactions IDs that are pseudoreactions with gpr.
    """

    def __init__(self, model: cobra.Model) -> None:

        # cobra version
        # basics
        self.name = model.id
        self.reac = model.reactions
        self.meta = model.metabolites
        self.gene = len(model.genes)
        # ends
        meta_ordedi = get_orphans_deadends_disconnected(model)
        self.orphans = meta_ordedi[0]
        self.deadends = meta_ordedi[1]
        self.disconnects = meta_ordedi[2]
        # balance
        mass_charge = get_mass_charge_unbalanced(model)
        self.mass_charge_unbalanced = [_ for _ in mass_charge[0] if _ in mass_charge[1]]
        self.mass_unbalanced = [
            _ for _ in mass_charge[0] if not _ in self.mass_charge_unbalanced
        ]
        self.charge_unbalanced = [
            _ for _ in mass_charge[1] if not _ in self.mass_charge_unbalanced
        ]
        # gpr
        self.pseudo = [_.id for _ in model.boundary] + test_biomass_presence(model)
        with_gpr = get_reac_with_gpr(model)
        self.normal_with_gpr = with_gpr[0]
        self.pseudo_with_gpr = with_gpr[1]

    def format_table(self, all_counts=True) -> pd.DataFrame:
        """Put the information of the report into a pandas DataFrame table.

        Args:
            - all_counts (bool, optional):
                Option to save the list of e.g. reactions
                as such or to convert them into counts when set to True.
                Defaults to True.

        Returns:
            pd.DataFrame:
                The data in table format
        """

        data = {
            "model": [self.name],
            "#reactions": [len(self.reac)],
            "#metabolites": [len(self.meta)],
            "#genes": [self.gene],
            "orphans": (
                [", ".join(self.orphans)] if not all_counts else [len(self.orphans)]
            ),
            "dead-ends": (
                [", ".join(self.deadends)] if not all_counts else [len(self.deadends)]
            ),
            "disconnects": (
                [", ".join(self.disconnects)]
                if not all_counts
                else [len(self.disconnects)]
            ),
            "mass and charge unbalanced": (
                [", ".join(self.mass_charge_unbalanced)]
                if not all_counts
                else [len(self.mass_charge_unbalanced)]
            ),
            "mass unbalanced": (
                [", ".join(self.mass_unbalanced)]
                if not all_counts
                else [len(self.mass_unbalanced)]
            ),
            "charge unbalanced": (
                [", ".join(self.charge_unbalanced)]
                if not all_counts
                else [len(self.charge_unbalanced)]
            ),
            "#normal reactions with gpr": (
                [", ".join(self.normal_with_gpr)]
                if not all_counts
                else [len(self.normal_with_gpr)]
            ),
            "#pseudoreactions with gpr": (
                [", ".join(self.pseudo_with_gpr)]
                if not all_counts
                else [len(self.pseudo_with_gpr)]
            ),
        }
        return pd.DataFrame(data)

    @implement
    def make_html():
        pass

    def visualise(self, color_palette: str = "YlGn") -> matplotlib.figure.Figure:
        """Visualise the basic information of the report.

        Args:
            - color_palette (str, optional):
                Colour palette to use for the plots.
                Defaults to 'YlGn'.

        Returns:
            matplotlib.figure.Figure:
                The visualisation as a single figure.
        """

        # basic settings
        # --------------

        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        # function to adjust position of labels
        def adjust_label_position(autotexts, min_dist):
            positions = np.array([text.get_position() for text in autotexts])
            for i in range(len(autotexts)):
                if autotexts[i].get_text() != "":
                    for j in range(i + 1, len(autotexts)):
                        if autotexts[i].get_text() != "":
                            dist = np.linalg.norm(positions[i] - positions[j])
                            if dist < min_dist:
                                shift = (min_dist - dist) / 2
                                delta = positions[j] - positions[i]
                                direction = delta / np.linalg.norm(delta)
                                positions[i] -= direction * shift
                                positions[j] += direction * shift

            for i, text in enumerate(autotexts):
                text.set_position(positions[i])

        # set up the figure
        fig = plt.figure()
        fig.suptitle(f"Basic information for model {self.name}", fontsize=16)
        grid = gspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.4)

        # 1: plot reacs, metabs and gene counts
        # -------------------------------------

        ax1 = fig.add_subplot(grid[0, 0])
        p = ax1.bar(
            ["reactions", "metabolites", "genes"],
            [len(self.reac), len(self.meta), self.gene],
            color=[cmap(0.25), cmap(0.5), cmap(0.8)],
            # edgecolor='black',
        )
        ax1.bar_label(p, [len(self.reac), len(self.meta), self.gene])
        ax1.set_ylabel("count")
        ax1.tick_params(axis="both", which="major", labelsize=9)
        ax1.set_title("A) Overview")
        ax1.set_ylim(top=ax1.get_ylim()[1] + ax1.get_ylim()[1] * 0.05)
        # ax1.set_xlabel('model entity')

        # 2: plot reacs with gpr
        # ----------------------

        local_grid = gspec.GridSpecFromSubplotSpec(
            1, 2, subplot_spec=grid[1, :], wspace=0.75
        )
        ax3 = plt.Subplot(fig, local_grid[0, 0])
        fig.add_subplot(ax3)
        ax4 = plt.Subplot(fig, local_grid[0, 1])
        fig.add_subplot(ax4)
        # ax3 = fig.add_subplot(grid[1,:])

        # plot reacs with gpr
        pie_data = [
            len(self.normal_with_gpr),
            (len(self.reac) - len(self.pseudo)) - len(self.normal_with_gpr),
            len(self.pseudo_with_gpr),
            len(self.pseudo) - len(self.pseudo_with_gpr),
        ]

        pie_label = ["normal/+", "normal/-", "pseudo/+", "pseudo/-"]

        def func(pct, allvals):
            absolute = int(np.round(pct / 100.0 * np.sum(allvals)))
            if absolute == 0:
                return ""
            return f"{pct:.1f}% ({absolute:d})"

        wedges, texts, autotexts = ax3.pie(
            pie_data,
            autopct=lambda pct: func(pct, pie_data),
            wedgeprops=dict(width=0.4),
            colors=[cmap(0.2), cmap(0.6), cmap(0.8), cmap(0.4)],
        )

        adjust_label_position(autotexts, 0.2)

        ax3.legend(
            wedges,
            pie_label,
            title="reaction +/- gpr",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
        )

        # plot reacs which are unbalanced
        pie_data = [
            len(self.mass_charge_unbalanced),
            len(self.charge_unbalanced),
            len(self.mass_unbalanced),
        ]

        pie_label = ["mass and charge", "charge only", "mass only"]

        def func(pct, allvals):
            absolute = int(np.round(pct / 100.0 * np.sum(allvals)))
            if absolute == 0:
                return ""
            return f"{pct:.1f}% ({absolute:d})"

        wedges, texts, autotexts = ax4.pie(
            pie_data,
            autopct=lambda pct: func(pct, pie_data),
            wedgeprops=dict(width=0.4),
            colors=[cmap(0.2), cmap(0.4), cmap(0.6)],
        )

        adjust_label_position(autotexts, 0.2)

        ax4.legend(
            wedges,
            pie_label,
            title="unbalanced",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
        )

        fig.text(0.5, 0.45, "C) Reactions", ha="center", va="center", fontsize=12)

        # 3: plot deadends, orhphans etc. for metabs
        # ------------------------------------------
        ax2 = fig.add_subplot(grid[0, 1])
        pie_data = [len(self.deadends), len(self.orphans), len(self.disconnects)]
        pie_data.append(len(self.meta) - sum(pie_data))

        pie_label = ["dead-ends", "orphans", "disconnects", "rest"]

        def func(pct, allvals):
            absolute = int(np.round(pct / 100.0 * np.sum(allvals)))
            if absolute == 0:
                return ""
            return f"{pct:.1f}% ({absolute:d})"

        wedges, texts, autotexts = ax2.pie(
            pie_data,
            autopct=lambda pct: func(pct, pie_data),
            explode=[0.1, 0.1, 0.1, 0],
            wedgeprops=dict(width=0.4),
            colors=[cmap(0.2), cmap(0.6), cmap(0.8), cmap(0.4)],
        )

        adjust_label_position(autotexts, 0.25)

        ax2.legend(
            wedges,
            pie_label,
            title="Metabolites",
            loc="upper left",
            bbox_to_anchor=(1, 0, 0.5, 1),
        )

        ax2.set_title("B) Metabolites")

        plt.setp(autotexts, size=8)

        return fig

    def save(self, dir: str, color_palette: str = "YlGn") -> None:
        """Save the report.

        Args:
            - dir (str):
                Directory to save the report to.
            - color_palette (str, optional):
                Colour palette of matplotlib to plot
                figures in. Defaults to 'YlGn'.
        """

        # save the statistics report
        self.format_table().to_csv(Path(dir, f"{self.name}_report.csv"), sep=";")

        # save the visualisation
        fig = self.visualise(color_palette)
        fig.savefig(Path(dir, f"{self.name}_visual.png"), bbox_inches="tight", dpi=400)

        # save the ids for unbalanced, gpr, ordedi
        balance = []
        for reaction in self.reac:
            if reaction.id in self.mass_charge_unbalanced:
                balance.append((reaction.id, "mass and charge unbalanced"))
            elif reaction.id in self.mass_unbalanced:
                balance.append((reaction.id, "mass unbalanced only"))
            elif reaction.id in self.charge_unbalanced:
                balance.append((reaction.id, "charge unbalanced only"))
            else:
                balance.append((reaction.id, "balanced"))
        pd.DataFrame(balance).to_csv(
            Path(dir, f"{self.name}_id_balance.csv"), sep=";", header=False
        )

        gpr = []
        for reaction in self.reac:
            if reaction.id in self.normal_with_gpr:
                gpr.append((reaction.id, "normal reaction with gpr"))
            elif reaction.id in self.pseudo_with_gpr:
                gpr.append((reaction.id, "pseudoreaction with gpr"))
            elif reaction.id in self.pseudo:
                gpr.append((reaction.id, "pseudoreaction without gpr"))
            else:
                gpr.append((reaction.id, "normal reaction without gpr"))
        pd.DataFrame(gpr).to_csv(
            Path(dir, f"{self.name}_id_gpr.csv"), sep=";", header=False
        )

        ordedi = []
        for metabolite in self.meta:
            if metabolite.id in self.orphans:
                ordedi.append((metabolite.id, "orphan"))
            elif metabolite.id in self.deadends:
                ordedi.append((metabolite.id, "deadend"))
            elif metabolite.id in self.disconnects:
                ordedi.append((metabolite.id, "disconnect"))
            else:
                ordedi.append((metabolite.id, "rest"))
        pd.DataFrame(ordedi).to_csv(
            Path(dir, f"{self.name}_id_ordedi.csv"), sep=";", header=False
        )


class MultiModelInfoReport(Report):

    def __init__(self) -> None:
        # super().__init__()
        self.table = pd.DataFrame(
            "model",
            "#reactions",
            "#metabolites",
            "#genes",
            "orphans",
            "dead-ends",
            "disconnects",
            "mass unbalanced",
            "charge unbalanced",
            "#reactions with gpr",
        )

    def add_single_report(self, report: ModelInfoReport) -> None:
        self.table = pd.concat([self.table, report], ignore_index=True)

    def __add__(self, other):
        self.table = pd.concat(self.table, other.table)

    @implement
    def visualise(self):
        pass

    @implement
    def save(self):
        pass


class GapFillerReport(Report):
    """Report for the gap-filling of the model.

    Attributes:
        - variety:
            The variety of the gap-filling method used.
        - statistics:
            List of different counts for reactions, genes and metabolites.
        - manual_curation:
            List of IDs for manual curation.
        - hide_zeros:
            Option to hide all zero values in the statistics. Defaults to False.
    """

    def __init__(
        self,
        variety: str,
        statistics: dict,
        manual_curation: dict,
        hide_zeros: bool = False,
        no_title: bool = False,
    ) -> None:
        self.variety = variety
        self.manual_curation = manual_curation
        self.hide_zeros = hide_zeros
        self.statistics = statistics
        self.no_title = no_title

    @property
    def statistics(self):
        """
        | Get or set the current statistics dictionary.
        | While setting the provided all zeros are removed from the dictionary
        | and if the values behind the keys 'unmappable' and 'missing (remaining)'
        | are the same only 'missing (remaining)' is kept.
        """
        return self._statistics

    @statistics.setter
    def statistics(self, stats_dict: dict):
        if self.hide_zeros:
            for k in stats_dict.keys():
                stats_dict[k] = {k: v for k, v in stats_dict[k].items() if v}

        self._statistics = stats_dict

    def visualise(self, color_palette: str = "YlGn") -> matplotlib.figure.Figure:
        """Visualise the basic information of the report.

        Args:
            - color_palette (str, optional):
                Colour palette to use for the plots.
                Defaults to 'YlGn'.

        Returns:
            matplotlib.figure.Figure:
                The visualisation as a single figure.
        """

        # basic settings
        # --------------

        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        fig = plt.figure(tight_layout=True)
        if not self.no_title:
            fig.suptitle(f"Statistics for Gapfilling via {self.variety}", fontsize=16)
        grid = gspec.GridSpec(2, 1, hspace=0.6)

        # plot statistics about genes
        # ---------------------------
        genes = list(self.statistics["genes"].keys())
        genes = genes[::-1]
        values = list(self.statistics["genes"].values())
        values = values[::-1]

        ax2 = fig.add_subplot(grid[0, 0])
        p = ax2.barh(
            genes, values, color=[cmap(0.2), cmap(0.4), cmap(0.6), cmap(0.8), cmap(1.0)]
        )
        ax2.bar_label(p, values)
        ax2.set_xlabel("count")
        ax2.tick_params(axis="x", which="major", labelsize=7, labelrotation=90)
        ax2.set_title("A) Genes")

        # plot statistics about reactions
        # -------------------------------
        reactions = list(self.statistics["reactions"].keys())
        reactions = reactions[::-1]
        values = list(self.statistics["reactions"].values())
        values = values[::-1]

        ax1 = fig.add_subplot(grid[1, 0])
        p = ax1.barh(
            reactions,
            values,
            color=[cmap(0.2), cmap(0.4), cmap(0.6), cmap(0.8), cmap(1.0)],
        )
        ax1.bar_label(p, values)
        ax1.set_xlabel("count")
        ax1.tick_params(axis="x", which="major", labelsize=7, labelrotation=90)
        ax1.set_title("B) Reactions")

        return fig

    def save(self, dir=str, color_palette: str = "YlGn") -> None:
        """Save the report.

        Args:
            - dir (str):
                Path to a directory to save the report to.
            - color_palette (str, optional):
                A colour gradient from the matplotlib library. If the name does not exist, uses the default. Defaults to `YlGn`.
        """
        dir_path = Path(
            dir,
            "GapFillerReport",
            f'{re.sub(r"[^a-zA-Z0-9]+", "_", self.variety).lower()}',
        )
        dir_path.mkdir(parents=True, exist_ok=True)

        # save the visualisation
        fig = self.visualise(color_palette)
        fig.savefig(
            Path(dir_path, "gapfill_statistics_visual.png"),
            bbox_inches="tight",
            dpi=400,
        )
        pd.DataFrame(self.statistics).to_csv(
            Path(dir_path, "gapfill_statistics.csv"), sep=";", header=True
        )

        # save the manual curation lists
        for category in self.manual_curation["reactions"]:
            if not self.manual_curation["reactions"][category].empty:
                pd.DataFrame(self.manual_curation["reactions"][category]).to_csv(
                    Path(dir_path, f"reactions_{category}.csv"),
                    sep=";",
                    header=True,
                    index=False,
                )
        for category in self.manual_curation["genes"]:
            if not self.manual_curation["genes"][category].empty:
                pd.DataFrame(self.manual_curation["genes"][category]).to_csv(
                    Path(dir_path, f"genes_{category}.csv"),
                    sep=";",
                    header=True,
                    index=False,
                )


class SBOTermReport(Report):
    """Report of the ABO terms of a model.

    Attributes:
        - name:
            Name (ID) of the model.
        - sbodata:
            Dictionary containing the SBO terms and the corresponding
            counts of annotations found in the model.
            Only includes SBO terms, that have at least 1 occurence in
            the model.
    """

    def __init__(self, model: libModel):
        """
        Args:
            - model (libModel):
                A model loaded with libSBML.
        """
        self.name = model.getId()
        self.sbodata = get_reactions_per_sbo(model)

    def visualise(self) -> matplotlib.figure.Figure:
        """Visualise the amount of SBO terms found in the model
        the report was created with.

        Returns:
            matplotlib.figure.Figure:
                The created graphic.
        """

        df = (
            pd.DataFrame(self.sbodata, index=[0])
            .T.reset_index()
            .rename({0: self.name, "index": "SBO-Term"}, axis=1)
        )
        df["SBO-Name"] = df["SBO-Term"].apply(search_sbo_label)
        ax = (
            df.drop("SBO-Term", axis=1)
            .sort_values(self.name)
            .set_index("SBO-Name")
            .plot.barh(width=0.8, figsize=(8, 10))
        )
        ax.set_ylabel("")
        ax.set_xlabel("number of reactions", fontsize=16)
        ax.legend(loc="lower right")
        fig = ax.get_figure()
        plt.tight_layout()

        return fig

    def save(self, dir: str):
        """Save the information inside

        Args:
            - dir (str):
                String to the output directory.
        """
        fig = self.visualise()
        fig.savefig(Path(dir, "sboterms.png"), dpi=400)


class MultiSBOTermReport:
    """A collection of SBO term reports.

    Attributes:
        - model_reports:
            List of :py:class:`~refinegems.classes.reports.SBOTermReports`.
    """

    def __init__(self, reports: Union[list[SBOTermReport] | SBOTermReport]):
        """
        Args:
            - reports (Union[list[SBOTermReport] | SBOTermReport]):
                Either a single or a list of
                :py:class:`~refinegems.classes.reports.SBOTermReports`

        Raises:
            - ValueError: Wrong input type
        """

        match reports:
            case list():
                self.model_reports = reports
            case SBOTermReport():
                self.model_reports = [reports]
            case _:
                raise ValueError(
                    "Parameter reports must be a list of SBOTermReports or a single SBOTermReport"
                )

    def add_report(self, report: SBOTermReport):
        """Add another :py:class:`~refinegems.classes.reports.SBOTermReports`
        to the report collection.

        Args:
            - report (SBOTermReport):
                The report to add.
        """
        self.model_reports.append(report)

    def visualise(
        self,
        rename: dict = Union[None, dict],
        color_palette: Union[str, list[str]] = "Paired",
        figsize: tuple = (10, 10),
    ) -> matplotlib.figure.Figure:
        """Visualise the amount of SBO terms in the models.

        Args:
            - rename (Union[None,dict], optional):
                Takes a dictioanry of model IDs and alternative names
                When set, uses the dictionary to rename the models.
                Defaults to None.
            - color_palette (Union[str,list[str]], optional):
                Color palette name or list of colours for the graphic.
                Defaults to 'Paired'.
            - figsize (tuple, optional):
                Site of the figure. Requires a tuple of two integers.
                Defaults to (10,10).

        Raises:
            - TypeError: Unkown type for color palette.

        Returns:
            matplotlib.figure.Figure:
                The generated graphic
        """

        # data formatting
        df = pd.DataFrame.from_dict({r.name: r.sbodata for r in self.model_reports})
        df = df.reset_index().rename({"index": "SBO-Term"}, axis=1)
        df["SBO-Name"] = df["SBO-Term"].apply(search_sbo_label)

        # set the colours
        match color_palette:
            case str():
                try:
                    cmap = matplotlib.colormaps[color_palette]
                except ValueError:
                    logging.WARN('Unknown color palette, setting it to "Paired"')
                    cmap = matplotlib.colormaps["Paired"]
                if isinstance(cmap, matplotlib.colors.ListedColormap):
                    cmap = cmap.colors[0 : len(self.model_reports)]
                else:
                    cmap = cmap(np.linspace(0, 1, len(self.model_reports)))
            case list():
                cmap = color_palette
            case _:
                mes = "Unkown type for color_palette"
                raise TypeError(mes)

        # prepare the data
        df.drop("SBO-Term", axis=1, inplace=True)
        df = df.set_index("SBO-Name")
        # sort by total SBO term count
        df["rowsum"] = df.sum(axis=1)
        df = df.sort_values(by="rowsum", axis=0)
        df.drop("rowsum", axis=1, inplace=True)

        # rename to custom names
        if rename is not None:
            df = df.rename(rename, axis=1)

        # make the figure
        ax = df.plot.barh(stacked=True, width=0.8, figsize=figsize, color=cmap)
        for patch in ax.patches:
            colour = patch.get_facecolor()
            patch.set_edgecolor(colour)
        ax.set_ylabel("")
        ax.set_xlabel("number of reactions", fontsize=16)
        ax.legend(loc="lower right")

        return ax.get_figure()

    def save(
        self,
        dir: str,
        rename: dict = Union[None, dict],
        color_palette: Union[str, list[str]] = "Paired",
        figsize: tuple = (10, 10),
    ):
        """Save the information of contained in the report.

        Args:
            - dir (str):
                String for the path of the outpt directory.
            - rename (Union[None,dict], optional):
                Takes a dictioanry of model IDs and alternative names
                When set, uses the dictionary to rename the models.
                Defaults to None.
            - color_palette (Union[str,list[str]], optional):
                Color palette name or list of colours for the graphic.
                Defaults to 'Paired'.
            - figsize (tuple, optional):
                Site of the figure. Requires a tuple of two integers.
                Defaults to (10,10).
        """
        fig = self.visualise(rename, color_palette, figsize)
        fig.save(Path(dir, "sboterms.png"), dpi=400)
