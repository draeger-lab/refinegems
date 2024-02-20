"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune, Famke Baeuerle, Gwendolyn O. Döbel'

################################################################################
# requirements
################################################################################

import cobra
import copy
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
import numpy as np
import pandas as pd
import seaborn as sns
import warnings

from importlib.resources import files
from itertools import chain
from pathlib import Path
from typing import Literal

from ..analysis.investigate import get_mass_charge_unbalanced, get_orphans_deadends_disconnected, get_num_reac_with_gpr

################################################################################
# variables
################################################################################

KEGG_GLOBAL_PATHWAY = {'01100': 'Metabolic pathways',
                       '01110': 'Biosynthesis of secondary metabolites',
                       '01120': 'Microbial metabolism in diverse environments'}

KEGG_OVERVIEW_PATHWAY = {'01200': 'Carbon metabolism',
                         '01210': '2-Oxocarboxylic acid metabolism',
                         '01212': 'Fatty acid metabolism',
                         '01230': 'Biosynthesis of amino acids',
                         '01232': 'Nucleotide metabolism',
                         '01250': 'Biosynthesis of nucleotide sugars',
                         '01240': 'Biosynthesis of cofactors',
                         '01220': 'Degradation of aromatic compounds'}

# @TODO # : Where to put this file // check connection (if it truly works)
KEGG_METABOLISM_PATHWAY = files('refinegems.data.pathway').joinpath('KEGG_pathway_metabolism.csv')
KEGG_METABOLISM_PATHWAY_DATE = "6. July 2023"

################################################################################
# classes
################################################################################

class Report():
    # @IDEA 
    # colour palette for visulisation
    # date/time of creation?
    pass


class SingleGrowthSimulationReport(Report):
    """Report for a single growth simulation, one media against one model.

    Attributes:
        model_name: Name of the model.
        medium_name: Name of the medium.
        growth_value: Simulated growth value.
        doubling_time: Simulated doubling time.
        additives: List of substances, that were added.
        no_exchange: List of substances that normally would be found in the media
            but have been removed, as they are not part of the exchange reactions
            of the model.
    """
    
    def __init__(self, model_name = None,
                 medium_name = None,
                 growth_value = None,
                 doubling_time = None,
                 additives = None,
                 no_exchange = None):
        self.model_name = model_name
        self.medium_name = medium_name
        self.growth_value = growth_value
        self.doubling_time = doubling_time
        self.additives = additives
        self.no_exchange = no_exchange


    def __str__(self):
        return (f"model: {self.model_name}\n"
                f"medium: {self.medium_name}\n"
                f"growth: {self.growth_value}\n"
                f"doubling time: {self.doubling_time}\n"
                f"additives: {self.additives}\n"
                f"no_exchange: {self.no_exchange}\n" 
        )
    

    def to_dict(self) -> dict:

        return {'model_name': self.model_name,
                 'medium_name': self.medium_name,
                 'growth_value': self.growth_value,
                 'doubling_time': self.doubling_time,
                 'additives': self.additives if len(self.additives) > 0 else None,
                 'no_exchange': self.no_exchange}


class GrowthSimulationReport(Report):
    """Report for the growth simulation analysis.

    Attributes:
        reports: List of the report for the single growth analysis.
        model: List of the model names.
        media: List of the media names.
    """

    def __init__(self, reports:list[SingleGrowthSimulationReport] = []):

        self.reports = reports
        self.models = set([_.model_name for _ in reports])
        self.media = set([_.medium_name for _ in reports])

    def __str__(self) -> str:
        
        return "\n\n".join(str(_) for _ in self.reports)

    def add_sim_results(self, new_rep: SingleGrowthSimulationReport):
        """Add a new single growth report to the reports list

        Args:
            new_rep (SingleGrowthSimulationReport): The new simulation report.
        """

        self.reports.append(new_rep)
        self.models.add(new_rep.model_name)
        self.media.add(new_rep.medium_name)

    def to_table(self) -> pd.DataFrame:
        """Return a table of the contents of the report.

        Returns:
            pd.DataFrame: The table containing the information in the report.
        """

        l = []
        for report in self.reports:
            l.append(report.to_dict())
        return pd.DataFrame(l)
    

    # @TODO
    # @NOTE: clean up for unrealistically high and minicules values to 0 - anyone a better idea?  
    def plot_growth(self, unit:Literal['h','dt']='dt', color_palette:str='YlGn') -> matplotlib.figure.Figure:
        """Visualise the contents of the report.

        Args:
            unit (Literal['h','dt'], optional): Set the unit to plot. 
                Can be doubling time in minutes ('dt') or growth rates in mmol/gDWh ('h'). 
                Defaults to 'dt'.
            color_palette (str, optional): A colour gradient from the matplotlib library.
                If the name does not exist, uses the default. 
                Defaults to 'YlGn'.

        Returns:
            matplotlib.figure.Figure: The plotted figure.
        """

        def plot_growth_bar(xdata:list[str], xlab:str, ydata:list[float], 
                            ylab:str, title:str, color_palette:str='YlGn') -> matplotlib.figure.Figure:
            """Helper function to plot the bar plot for the growth visualisation.

            Args:
                xdata (list[str]): List of the x-axis data (medium or model names).
                xlab (str): The x-axis label.
                ydata (list[float]): List of thr y-axis data (the values).
                ylab (str): The y-axis label.
                title (str): The title of the plot.
                color_palette (str, optional): A colour gradient from the matplotlib library.
                    If the name does not exist, uses the default. 
                    Defaults to 'YlGn'.

            Returns:
                matplotlib.figure.Figure: The plotted figure.
            """
            
            # create colour gradient
            try:
                cmap = matplotlib.colormaps[color_palette]
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.colormaps['YlGn']

            # set up the figure
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])

            # clean-up data
            # @TODO / @NOTE
            ydata = [_ if _ > 0.0 else 0.0 for _ in ydata]
            ydata = [_ if _ < 1000.0 else 0.0 for _ in ydata]

            # construct the plot
            max_ydata = max(ydata)
            if max_ydata <= 0:
                warnings.warn('Model is not able to grow sensible on any medium. Returning empty figure.')
                return fig
            cont = ax.bar(list(xdata), ydata, color=cmap([_/max_ydata for _ in ydata]))

            # set labels and others
            ax.bar_label(cont, fmt='%.2f', color='black', padding=1.0)
            ax.set_ylabel(ylab, labelpad=12)
            ax.set_xlabel(xlab, labelpad=12)
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
            ax.set_title(title)

            return fig


        def plot_growth_heatmap(data: pd.DataFrame, color_palette:str='YlGn') -> matplotlib.figure.Figure:
            """Helper function to plot the heatmap for the growth visualisation.

            Args:
                data (pd.DataFrame): The table containing the data to be plotted.
                    Needs to have the columns 'medium', 'model' and one for the growth values.
                color_palette (str, optional): A colour gradient from the matplotlib library.
                    If the name does not exist, uses the default. 
                    Defaults to 'YlGn'.

            Returns:
                matplotlib.figure.Figure: The plotted figure.
            """

            # clean up + transform data
            growth=data.set_index(['medium', 'model']).sort_index().T.stack()
            growth.columns.name=None
            growth.index.names = (None,None)
            growth.index.name=None
            growth.index = growth.index.get_level_values(1)
            growth[growth > 500] = 0
            growth[growth < 0] = 0
            growth.replace([np.inf, -np.inf], 0, inplace=True)
            over_growth = growth.max().max() + 6
            growth.replace(np.nan, over_growth, inplace=True)
            under_growth = growth.min().min() - 5
            vmin= under_growth if under_growth > 1e-5 else 1e-5 #Use same threshhold as in find_missing_essential in growth
            vmax=over_growth - 1
            annot = growth.copy()
            annot = annot.round().astype(int)
            annot[annot < 1e-5] = ''
            annot.replace(over_growth.round().astype(int), 'No data', inplace=True)

            # setting the colours 
            try:
                cmap = matplotlib.colormaps[color_palette]
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.colormaps['YlGn']
            cmap.set_under('black') # too low / no growth
            cmap.set_over('white') # no data

            # plot the heatmap
            fig, ax = plt.subplots(figsize=(10,8))
            sns.heatmap(growth.T, 
                        annot=annot.T, 
                        annot_kws={"fontsize":15},
                        vmin=vmin, 
                        vmax=vmax,
                        cmap=cmap, 
                        linewidth=.5, 
                        cbar_kws = {'orientation':'vertical', 'label':'doubling time [min]', 'extend': 'min', 'extendrect':True},
                        ax=ax,
                        fmt=''
                        )
            rotation = 40 if len(growth.index) > 3 else 0
            plt.tick_params(rotation=0, bottom=False, top=False, left=False, right=False)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha="right")

            return fig
        

        # match the unit 
        match unit:
            case 'h':
                unit_text = r'growth rate $[\frac{mmol}{gDWh}]$'
            case 'dt':
                unit_text = 'doubling time [min]'
            case _:
                raise ValueError(f'Unknown input for unit: {unit}')

        # multiple models vs one medium
        if len(self.models) > 1 and len(self.media) == 1:
            
            # collect data 
            xdata = [_.model_name for _ in self.reports]
            xlab = 'models'
            ydata = [_.growth_value if unit=='h' else _.doubling_time for _ in self.reports]
            ylab = unit_text
            title = f'Growth simulation on {next(iter(self.media))} for different models'

            # plot
            return plot_growth_bar(xdata,xlab,ydata,ylab,title, color_palette)

        # one model vs mutiple media
        elif len(self.models) == 1 and len(self.media) > 1:

            # collect data
            xdata = [_.medium_name for _ in self.reports]
            xlab = 'media'
            ydata = [_.growth_value if unit=='h' else _.doubling_time for _ in self.reports]
            ylab = unit_text
            title = f'Growth simulation for {next(iter(self.models))} on different media'

            # plot
            return plot_growth_bar(xdata,xlab,ydata,ylab,title, color_palette)
        
        # one medium, one model case - just to make it usable for all inputs
        elif len(self.models) == 1 and len(self.media) == 1:

            # collect data
            xdata = [_.medium_name for _ in self.reports]
            xlab = 'medium'
            ydata = [_.growth_value if unit=='h' else _.doubling_time for _ in self.reports]
            ylab = unit_text
            title = f'Growth simulation for {next(iter(self.models))} on medium {xdata[0]}'

            # plot
            return plot_growth_bar(xdata,xlab,ydata,ylab,title, color_palette)

        # multiple vs multiple
        elif len(self.models) > 1 and len(self.media) > 1:
            
            data = pd.DataFrame({'model':[_.model_name for _ in self.reports], 'medium':[_.medium_name for _ in self.reports], 'growth':[_.growth_value if unit=='h' else _.doubling_time for _ in self.reports]})

            return plot_growth_heatmap(data, color_palette)
        
        # problematic case
        else:
            raise IndexError('Can only plot growth if at least one model and one medium are present.')

        
    # @TEST    
    # @EXTEND : more options for saving the report e.g. html or pdf
    def save(self, to:str, how:Literal['dir']='dir', check_overwrite:bool=True, color_palette:str='YlGn'):
        """Save the report. Current options include:

        - 'dir': save the report to a directory, including a txt and two graphics
        - .... see future updates .....

        Args:
            to (str): Path to a directory to save the report to.
            how (Literal['dir'], optional): How to save the report. 
                For options see functions description. 
                Defaults to 'dir'.
            check_overwrite (bool, optional): Flag to choose to check for existing directory/files of same name 
                or just to overwrite them. Defaults to True.
            color_palette (str, optional): A colour gradient from the matplotlib library.
                If the name does not exist, uses the default. 
                Defaults to 'YlGn'.

        Raises:
            ValueError: If the parameter 'how' is given something unexpected.
        """

        match how:
            # save to a new directory
            case 'dir':
                # make sure given directory path ends with '/'
                    if not to.endswith('/'):
                        to = to + '/'
                    # create directory to save report to
                    dir_path = to + 'GrowthSimReport/'
                    Path(dir_path).mkdir(parents=True, exist_ok=check_overwrite)
                    # save the report
                    with open(dir_path + 'report.txt', 'w') as f:
                        f.write(str(self))
                    # save visualisation for doubling time
                    fig_dt = self.plot_growth(color_palette=color_palette)
                    fig_dt.savefig(dir_path + 'report_vis_dt.png', bbox_inches='tight')
                    # save visualisation for growth rate 
                    fig_dt = self.plot_growth(unit='h', color_palette=color_palette)
                    fig_dt.savefig(dir_path + 'report_vis_h.png', bbox_inches='tight')

            case _:
                raise ValueError(f'Unknow input for parameter "how": {how}.\n Cannot save report. Abort.')
            

class KEGGPathwayAnalysisReport(Report):
    """Report for the KEGG pathway analysis.

    Attributes:
        total_reac: An integer for the total number of reactions in the model.
        kegg_count: An integer as a counter for the KEGG pathway annotations.
        kegg_global: Dictionary of global KEGG IDs and their counts.
        kegg_over: Dictionary of overvire KEGG IDs and their counts.
        kegg_rest: Dictionary of the remaining KEGG IDs and their counts.
    """
    
    def __init__(self, 
                 total_reac=None, kegg_count=None,
                 kegg_global=None, kegg_over=None, kegg_rest=None) -> None:
        
        # super().__init__()
        # general counts
        self.total_reac = total_reac
        self.kegg_count = kegg_count

        # kegg pathways
        self.kegg_global = kegg_global
        self.kegg_over = kegg_over
        self.kegg_paths = kegg_rest


    def visualise_kegg_counts(self, colors:list[str]=['lightgreen','darkgreen']) -> plt.figure:
        """Visualise the amounts of reaction with and without
        KEGG pathway annotation.

        Args:
            colors (list[str], optional): List of two colours used for the plotting.
                If wrong number or non-matplotlib colours are given, sets its to the default.
                Defaults to 'lightgreen' and 'darkgreen'.

        Returns:
            plt.figure: The resulting plot.
        """

        # validate colors
        if len(colors) != 2 or not matplotlib.colors.is_color_like(colors[0]) or not matplotlib.colors.is_color_like(colors[1]):
            warnings.warn('Unknown colors or false amount of colors for pie chart. Resume using default values.')
            colors = ['lightgreen','lightskyblue']
        
        # generate the plot
        explode = (0.0,0.1)
        fig, ax = plt.subplots()
        values = [self.kegg_count, self.total_reac-self.kegg_count]
        labels = ['yes','no']
        ax.pie(values,
               autopct=lambda pct: "{:.1f}%\n({:.0f})".format(pct, (pct/100)*sum(values)),
               colors=colors,
               explode=explode, shadow = True, startangle=90)
        ax.legend(labels, title='KEGG\npathway')

        return fig


    def visualise_kegg_pathway(self, plot_type:Literal['global','overview','high','existing']='global', 
                               label:Literal['id','name']='id',
                               color_palette:str='YlGn') -> plt.figure:
        """Visualise the KEGG pathway identifiers present.

        Depending on the :plot_type:, different levels of pathway identifiers
        are plotted:

        - global: check and plot only the global pathway identifiers
        - overview: check and plot only the overview pathway identifiers
        - high: check and plot all identifiers grouped by their high level pathway identifiers. This option uses label=name, independedly of the input
        - all: check and plot all identifiers

        Args:
            plot_type (Literal["global","overview","high","existing"], optional): Type of plot, explaination see above. Defaults to 'global'.
            label (Literal["id","name"], optional): Type of the label. If 'id', uses the KEGG pathway IDs,
                if 'name', uses the pathway names. Defaults to 'id'.
            color_palette (str, optional): A colour gradient from the matplotlib library.
                If the name does not exist, uses the default. 
                Defaults to 'YlGn'.

        Returns:
            plt.figure: The plotted visualisation.
        """

        # get data and KEGG pathway label mapping
        # for the given plot type
        match plot_type:
            case 'global':
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = 'global identifiers'
            case 'overview':
                data = self.kegg_over
                label_map = KEGG_OVERVIEW_PATHWAY
                title_type = 'overview identifiers'
            case 'high':
                label = 'name'
                data = self.kegg_paths
                label_map = pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str).set_index('id')[['group']].to_dict()['group']
                title_type = 'grouped identifiers'
            case 'existing':
                data = self.kegg_paths
                label_map = pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str).set_index('id')[['specific']].to_dict()['specific']
                title_type = 'identifiers in model'
            case _:
                warnings.warn(F'Unknown option for plot_type, choosing "global" istead: {plot_type}')
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = 'global identifiers'

        # create the label
        match label:
            case 'id':
                for k in label_map:
                    if k not in data and plot_type != 'existing':
                        data[k] = 0
            case 'name':
                old_data = data
                data = {}
                for k in label_map:
                    if k not in old_data:
                        if plot_type != 'existing' and label_map[k] not in data:
                            data[label_map[k]] = 0
                    else:
                        if label_map[k] not in data:
                            data[label_map[k]] = old_data[k]
                        else:
                            data[label_map[k]] += old_data[k]

            case _:
                warnings.warn(F'Unknown input for label: {label}. Using "id" instead.')
                for k in label_map:
                    if k not in data:
                        data[k] = 0

        data = pd.DataFrame(data.items(), columns=['label', 'counts']).sort_values('counts')
        cdata = data.counts.values
        ldata = data.label.values

        # create the graph
        # ----------------
        # setting the colours 
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps['YlGn']
        # set up the figure
        fig = plt.figure()
        if 'existing' == plot_type:
            ax = fig.add_axes([0,0,5,6])
            # construct the plot
            cont = ax.barh(ldata, cdata, color=cmap([x/max(cdata) for x in cdata]),
                   label=ldata)
            ax.bar_label(cont, fmt='%d', color='black', padding=1.0)
            ax.set_ylabel('KEGG pathway')
            ax.set_xlabel('Number of reaction annotations')
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
            ax.set_title(F'Pathway analysis with KEGG: {title_type}')
        else:

            ax = fig.add_axes([0,0,1,1])
            # construct the plot
            cont = ax.barh(ldata, cdata, color=cmap([x/max(cdata) for x in cdata]),
                   label=ldata)
            ax.bar_label(cont, fmt='%d', color='black', padding=1.0)
            ax.set_ylabel('KEGG pathway', labelpad=12)
            ax.set_xlabel('Number of reaction annotations', labelpad=12)
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
            ax.set_title(F'Pathway analysis with KEGG: {title_type}')

        return fig


    def save(self, dir:str) -> None:
        """Save the content of the report as plots.

        Args:
            dir (str): Path to a directory to save the output directory with all the plot in.
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}pathway-analysis/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}pathway-analysis/"}')
        except FileExistsError:
            print('Given directory already has required structure.')

        # create and save plots
        # a) for the counts
        if self.total_reac and self.kegg_count:
            count_fig = self.visualise_kegg_counts()
            count_fig.savefig(F'{dir}pathway-analysis/kegg_anno_counts.png', bbox_inches='tight')
        # b) for the actual pathways
        # 1.) global KEGG IDs
        if self.kegg_global:
            # with id
            fig = self.visualise_kegg_pathway(plot_type='global', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_global_id.png', bbox_inches='tight')
            # with name
            fig = self.visualise_kegg_pathway(plot_type='global', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_global_name.png', bbox_inches='tight')
        # 2.) Overview KEGG IDs
        if self.kegg_over:
            # with id
            fig = self.visualise_kegg_pathway(plot_type='overview', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_overview_id.png', bbox_inches='tight')
            # with name
            fig = self.visualise_kegg_pathway(plot_type='overview', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_overview_name.png', bbox_inches='tight')
        # 3.) rest
        if self.kegg_paths:
            # grouped by high-level terms
            fig = self.visualise_kegg_pathway(plot_type='high', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_high.png', bbox_inches='tight')
            # all with id
            fig = self.visualise_kegg_pathway(plot_type='existing', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_existing_id.png', bbox_inches='tight')
            # all with name
            fig = self.visualise_kegg_pathway(plot_type='existing', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_existing_name.png', bbox_inches='tight')


class AuxotrophySimulationReport(Report):
    """Report for the auxotrophy simulation.

    Attributes:
        simulation_results: The data of the simulation.
    """
    
    def __init__(self, results) -> None:
        # super().__init__()
        self.simulation_results = results

    
    # @TEST
    # auxotrophy sim visualisation
    def visualise_auxotrophies(self, color_palette:str='YlGn', save:None|str=None) -> None|matplotlib.figure.Figure:
        """Visualise and/or save the results of the :py:func:`test_auxotrophies` function.

        Args:
            res (pd.DataFrame): The output of  :py:func:`test_auxotrophies`.
            color_palette (str, optional): A name of a seaborn gradient color palette. 
                In case name is unknown, takes the default. Defaults to 'YlGn'.
            save (None | str, optional): Path to a directory, if the output shall be saved. Defaults to None (returns the figure).

        Returns:
            None|matplotlib.figure.Figure: Either saves the figure and a table of the results or returns the plotted figure.
        """
        
        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps['YlGn']

        # set up the figure
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])

        # create heatmap
        sns.heatmap(self.simulation_results, ax=ax, cmap=cmap, cbar_kws={'label': 'flux'}, annot = True, fmt='.2f')

        # add labels
        ax.set_ylabel('amino acid', labelpad=12)
        ax.set_xlabel('medium', labelpad=12)
        ax.set_title('Fluxes for auxotrophy tests')

        # save or return
        if save:
            # make sure given directory path ends with '/'
            if not save.endswith('/'):
                save = save + '/'
            # save the visualisation of the growth rates
            fig.savefig(F'{save}auxotrophies_vis.png', bbox_inches='tight')
        
        else:
            return fig 
        
        
    def save(self, dir:str, color_palette:str='YlGn'):
        """Save the report to a given dictionary.

        Args:
            dir (str): Path to a dictionary.
            color_palette (str, optional): Name of a matplotlib colour palette. Defaults to 'YnGr'.
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'
        
        # save the visualisation of the growth rates
        self.visualise_auxotrophies(color_palette, save=dir)

        # save the growth rates as tabular information
        self.simulation_results.to_csv(F'{dir}auxotrophies_table.tsv', sep='\t', index=True)


class SourceTestReport(Report):
    """Report for the source test (:py:func:rg.growth.test_growth_with_source).

    Attributes:
        results: A pd.DataFrame with the results (substances and growth values).
        element: The element the test was performed for.
        model_name: The name of the model, that was tested.
    """
    
    def __init__(self, results:pd.DataFrame=None, element:str=None, model_name:str=None):
        # super().__init__()
        self.results = results
        self.element = element
        self.model_name = model_name


    def visualise(self, width:int=12, color_palette:str='YlGn') -> tuple[matplotlib.figure.Figure, pd.DataFrame]:
        """Visuale the results of the source test as a heatmap

        Args:
            width (int, optional): number of columns to display for the heatmap. 
                Number of row is calculated accordingly to fit all values.
                Defaults to 12.
            color_palette (str, optional): Color palette (gradient) for the plot. 
                Defaults to 'YlGn'.

        Returns:
            tuple(matplotlib.Figure, pd.DataFrame): The heatmap and the legend explaining the heatmap.
        """
        
        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps['YlGn']

        cmap.set_under('black') # too low / no growth
        cmap.set_over('white') # no data

        # get size of heatmap
        height = math.ceil(len(self.results)/width)
        total_cells = width * height

        # create table for plotting
        data_to_plot = copy.deepcopy(self.results)
        if len(self.results) < total_cells:
            temp = pd.DataFrame.from_records([['empty',None]]*(total_cells-len(self.results)),columns=['substance','growth value'])
            data_to_plot = pd.concat([data_to_plot, temp], ignore_index=True)
        data_to_plot['row'] = list(range(1,height+1))*width
        data_to_plot['column'] = list(chain.from_iterable([[x]*height for x in range(1,width+1)]))

        # remove unplottable entries 
        data_to_plot['growth value'].replace([np.inf, -np.inf], 0, inplace=True)
        over_growth = data_to_plot['growth value'].max() + 0.1 * data_to_plot['growth value'].max()
        data_to_plot['growth value'].replace(np.nan, over_growth, inplace=True)
        vmin= 1e-5 #Use same threshhold as in find_missing_essential in growth
        vmax=over_growth - 0.05 * data_to_plot['growth value'].max()

        # set annotations
        annot = data_to_plot.copy()
        annot['growth value'] = annot['growth value'].round(2)
        annot['growth value'] = annot['growth value'].apply(lambda x: '' if x < 1e-5 else (x if x < vmax else 'X'))

        detected_growth = len(annot[(annot['growth value'] != '') & (annot['growth value'] != 'X')])

        annot = annot.pivot(index='row',columns='column', values='growth value')
        legend = data_to_plot.pivot(index='row',columns='column', values='substance')

        # plot
        ax = sns.heatmap(data_to_plot.pivot(index='row',columns='column', values='growth value'),
                        linewidth=.5, cmap=cmap,
                        vmin=vmin, vmax=vmax,
                        annot=annot, fmt='',
                        cbar_kws={'label': r'growth rate $[\frac{mmol}{gDWh}]$'}
                        )

        ax.set(xlabel='column', ylabel='row')
        ax.set_title(f'Growth detected with {detected_growth} sources', fontsize=10)
        plt.suptitle(f'{self.element}-source growth simulation on model {self.model_name}')

        return (ax.get_figure(), legend)


    def save(self, dir:str, width:int=12, color_palette:str='YlGn') -> None:
        """Save the results of the source test.

        Args:
            dir (str): Path to a directory to save the results to.
            width (int, optional): Number of columns for the heatmap. 
                Defaults to 12.
            color_palette (str, optional):Color palette (gradient) for the plot. 
                Defaults to 'YlGn'.
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'
        
        # save the list 
        self.results.to_csv(dir+'source_test_results.csv', sep=';', header=True, index=False)

        # save the visualisation
        fig,leg = self.visualise(width=width, color_palette=color_palette)
        fig.savefig(dir+'source_test_hm.png', bbox_inches='tight', dpi=400)
        leg.to_csv(dir+'source_test_hm_legend.csv', sep=';', header=True, index=True)


class CorePanAnalysisReport(Report):
    """Report for the core-pan analysis. 

    Summarises the information and provides functions 
    for visualisation.

    Attributes:
        model: The model the report is based on.
        core_reac: List of reactions considered "core".
        pan_reac: List of reactions considered "pan".
        novel_reac: List of reactions considered "novel".
    """

    def __init__(self, model: cobra.Model,
                 core_reac:list[str]=None, pan_reac:list[str]=None, novel_reac:list[str]=None):

        # super().__init__()
        # general attributes
        self.model = model
        # reaction attributes
        self.core_reac = core_reac
        self.pan_reac = pan_reac
        self.novel_reac = novel_reac
        # ...

    def get_reac_counts(self):
        """Return a dictionary of the counts of the reactions types (core, pan, novel).
        """

        counts = {}
        if self.core_reac:
            counts['core'] = len(self.core_reac)
        else:
            counts['core'] = 0
        if self.pan_reac:
            counts['pan'] = len(self.pan_reac)
        else:
            counts['pan'] = 0
        if self.novel_reac:
            counts['novel'] = len(self.novel_reac)
        else:
            counts['novel'] = 0

        return counts


    #@TODO
    def isValid(self,check='reaction-count') -> bool:
        """Check if a certain part of the analysis is valid.

        Currently possible checks:
            reaction-count : check if the number of reactions in the model
                             equal the sum of the novel, pan and core reactions

        @TODO
            implements more checks

        Args:
            check (str, optional): Describes which part to check. Options are listed above.
                Defaults to 'reaction-count'.

        Raises:
            ValueError: Unknown string for parameter check. 

        Returns:
            bool: Result of the check.
        """

        match check:
            case 'reaction-count':
                pc_total = sum(self.get_reac_counts().values())
                diff = len(self.model.reactions) - pc_total
                if diff != 0:
                    return False
                else:
                    return True
            case _:
                raise ValueError('Unknown string for parameter check: ', check)


    def visualise_reactions(self) -> matplotlib.figure:
        """Visualise the results of the pan-core analysis for the reactions as a donut chart.

        Returns:
            matplotlib.figure: The plot.
        """
        

        # check the counts
        # ----------------
        counts = self.get_reac_counts()

        if self.isValid('reaction-count'):
            vis_data = list(counts.values())
            vis_label = list(counts.keys())
            vis_color = ['lightgreen','lightskyblue','lightcoral']

        else:
            warnings.warn('Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected.')
            vis_data = list(counts.values()).append(len(self.model.reactions) - sum(counts.values))
            vis_label = list(counts.keys()).append('discrepancies')
            vis_color = ['lightgreen','lightskyblue','lightcoral', 'gold']

        # plot a donut chart
        # ------------------
        fig, ax = plt.subplots()
        wedges, texts, autotexts = ax.pie(vis_data,
                                          autopct=lambda pct: "{:.1f}%\n({:.0f})".format(pct, (pct/100)*sum(vis_data)),
                                          pctdistance=.8, labeldistance=1.25,
                                          colors=vis_color,
                                          radius = 1, wedgeprops=dict(width=0.4, edgecolor='w'))

        ax.legend(wedges, vis_label,
                  title="Classification",
                  loc="center left",
                  bbox_to_anchor=(1, 0, 0.5, 1))

        ax.set_title(F"Results of core-pan analysis\nof the reactions of model {self.model.id}")

        return fig


    #@TODO
    def save(self, dir:str):
        """Save the results inside a PanCoreAnalysisReport object.

        The function creates a new folder 'pan-core-analysis'
        inside the given directory and creates the following documents:

        - table_reactions.tsv : reactions ID mapped to their labels
        - visualise_reactions : donut chart of the values above

        Args:
            dir (str): Path to a directory to save the output to.
        """
        
        # ..........................................................
        #@TODO
        #    - an easily human readable overview file?
        #    - smth about metabolites?
        # ..........................................................

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}pan-core-analysis/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}pan-core-analysis/"}')
        except FileExistsError:
            print('Given directory already has required structure.')

        # save the reactions visualisation
        reac_vis = self.visualise_reactions()
        reac_vis.savefig(F'{dir}pan-core-analysis/visualise_reactions.png', dpi=reac_vis.dpi)

        # save table of reactions mapped to characterisation
        if not self.isValid(check='reaction-count'):
            warnings.warn('Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected. Only labbeld reactions will be written to table.')
        reac_tab = pd.DataFrame({'reaction_id': self.core_reac + self.pan_reac + self.novel_reac,
                                'pan-core': (['core']*len(self.core_reac)) + (['pan']*len(self.pan_reac)) + (['core']*len(self.novel_reac))})
        reac_tab.to_csv(F'{dir}pan-core-analysis/table_reactions.tsv', sep='\t', index=False)


# @TODO
class ModelInfoReport(Report):
    """Report about the basic information of a given model.

    Note: currently requires the input model to be a COBRApy model object.

    Attributes:
        name: A string for the name of the model.
        reac: An int that describes the number of reactions in the model.
        meta: An int that describes the number of metabolites in the model.
        gene: An int that describes the numver of genes in the model.
        orphans: List of metabolite IDs that are considered orphans.
        deadends: List of metabolite IDs that are considered dead-ends.
        disconnects: List of metabolites that are disconnected in the model.
        mass_unbalanced: List of reaction IDs that are unbalanced regarding their mass.
        charge_unbalanced: List of reactions IDs that are unbalanced regarding their charges.
        with_gpr: Integer describing the number of reactions that have a gene production rule.
    """
    
    def __init__(self, model) -> None:
        
        # cobra version
        # basics
        self.name = model.id
        self.reac = len(model.reactions)
        self.meta = len(model.metabolites)
        self.gene = len(model.genes)
        # ends
        meta_ordedi = get_orphans_deadends_disconnected(model)
        self.orphans = meta_ordedi[0]
        self.deadends = meta_ordedi[1]
        self.disconnects = meta_ordedi[2]
        # balance
        mass_charge = get_mass_charge_unbalanced(model)
        self.mass_unbalanced = mass_charge[0]
        self.charge_unbalanced = mass_charge[1]
        # gpr
        self.with_gpr = get_num_reac_with_gpr(model)

    def format_table(self, all_counts=True) -> pd.DataFrame:
        """Put the information of the report into a pandas DataFrame table.

        Args:
            all_counts (bool, optional): Option to save the list of e.g. reactions
                as such or to convert them into counts when set to True. 
                Defaults to True.

        Returns:
            pd.DataFrame: The data in table format
        """

        data = {'model': [self.name],
                '#reactions': [self.reac],
                '#metabolites': [self.meta],
                '#genes': [self.gene],
                'orphans': [', '.join(self.orphans)] if not all_counts else [len(self.orphans)],
                'dead-ends': [', '.join(self.deadends)] if not all_counts else [len(self.deadends)],
                'disconnects': [', '.join(self.disconnects)] if not all_counts else [len(self.disconnects)],
                'mass unbalanced': [', '.join(self.mass_unbalanced)] if not all_counts else [len(self.mass_unbalanced)],
                'charge unbalanced': [', '.join(self.charge_unbalanced)] if not all_counts else [len(self.charge_unbalanced)],
                '#reactions with gpr': [self.with_gpr]
                } 
        return pd.DataFrame(data)

    # @TODO
    def make_html():
        pass

    def visualise(self, color_palette:str='YlGn') -> matplotlib.figure.Figure:
        """Visualise the basic information of the report.

        Args:
            color_palette (str, optional): Colour palette to use for the plots. 
                Defaults to 'YlGn'.

        Returns:
            matplotlib.figure.Figure: The visualisation as a single figure.
        """

        # basic settings
        # --------------

        # create colour gradient
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps['YlGn']

        # set up the figure
        fig = plt.figure()
        fig.suptitle(f'Basic information for model {self.name}', fontsize=16)
        grid = gspec.GridSpec(2,2, height_ratios=[2,1], hspace=0.4)

        # 1: plot reacs, metabs and gene counts
        # -------------------------------------

        ax1 = fig.add_subplot(grid[0,0])
        p = ax1.bar(['reactions','metabolites','genes'],
            [self.reac,self.meta,self.gene],
            color=[cmap(0.25),cmap(0.5),cmap(0.8)]
            # edgecolor='black',
            )
        ax1.bar_label(p, [self.reac,self.meta,self.gene])
        ax1.set_ylabel('count')
        ax1.tick_params(axis='both', which='major', labelsize=9)
        ax1.set_title('A) Overview')
        ax1.set_ylim(top=ax1.get_ylim()[1] + ax1.get_ylim()[1]*0.05)
        # ax1.set_xlabel('model entity')

        # 2: plot reacs with gpr
        # ----------------------

        local_grid = gspec.GridSpecFromSubplotSpec(2,1, subplot_spec=grid[1,:])
        ax3 = plt.Subplot(fig, local_grid[0,0])
        fig.add_subplot(ax3)
        ax4 = plt.Subplot(fig, local_grid[1,0])
        fig.add_subplot(ax4)
        # ax3 = fig.add_subplot(grid[1,:])

        # plot reacs with gpr
        stacked_bars = {'with gpr': np.array([self.with_gpr]),
                        'no gpr': ([self.reac - self.with_gpr])}
        bottom = np.zeros(1)

        c = 0.3

        for label,count in stacked_bars.items():
            p = ax3.barh(['reactions'],count,
                        label=label, left=0.0,
                        color=[cmap(c)])
            ax3.bar_label(p, count, rotation=270)
            bottom += count
            c += 0.5

        ax3.set_title('C) Reactions')
        ax3.set_ylabel('gpr')
        ax3.tick_params(left = False,labelleft = False ,
                                labelbottom = False, bottom = False)
        ax3.legend(bbox_to_anchor=(0.75, 0, 0.5, 1.05), loc="center right")

        # plot reacs which are unbalanced
        mass_and_charge = [_ for _ in self.mass_unbalanced if _ in self.charge_unbalanced]
        only_mass = [_ for _ in self.mass_unbalanced if not _ in mass_and_charge]
        only_charge = [_ for _ in self.charge_unbalanced if not _ in mass_and_charge]

        stacked_bars = {'mass and charge': np.array([len(mass_and_charge)]),
                        'mass only': np.array([len(only_mass)]),
                        'charge only': np.array([len(only_charge)])}
        bottom = np.zeros(1)

        c = 0.3

        for label,count in stacked_bars.items():
            if count > 0:
                p = ax4.barh(['reactions'],count,
                            label=label, left=0.0,
                            color=[cmap(c)])
                ax4.bar_label(p, count, rotation=270)
                bottom += count
                c += 0.3

        ax4.set_xlabel('count')
        ax4.set_ylabel('unbal.')
        ax4.tick_params(left = False,labelleft = False ,
                                labelbottom = False, bottom = False)
        ax4.legend(title='unbalanced', 
                bbox_to_anchor=(0.875, 0,0.5, 0.25), loc="center right")



        # 3: plot deadends, orhphans etc. for metabs
        # ------------------------------------------
        ax2 = fig.add_subplot(grid[0,1])
        pie_data = [len(self.deadends), len(self.orphans), len(self.disconnects)]
        pie_data.append(self.meta - sum(pie_data))

        pie_label = ['dead-ends','orphans','disconnects','rest']

        def func(pct, allvals):
            absolute = int(np.round(pct/100.*np.sum(allvals)))
            return f"{pct:.1f}% ({absolute:d})"

        wedges, texts, autotexts  = ax2.pie(pie_data, autopct=lambda pct: func(pct, pie_data), 
            explode = [0.1, 0.1, 0.1, 0], wedgeprops=dict(width=0.4), 
            colors=[cmap(0.2),cmap(0.6), cmap(0.8), cmap(0.4)])

        ax2.legend(wedges, pie_label,
                title="Metabolites",
                loc="upper left",
                bbox_to_anchor=(1, 0, 0.5, 1)
                )

        ax2.set_title('B) Metabolites')

        plt.setp(autotexts, size=8)

        return fig


    # @TODO: match case for different output?
    #        -> only needed when html options available
    def save(self, dir:str, color_palette:str='YlGn') -> None:
        """Save the report.

        Args:
            dir (str): Directory to save the report to.
            color_palette (str, optional): Colour palette of matplotlib to plot
                figures in. Defaults to 'YlGn'.
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # save the statistics report
        self.format_table().to_csv(f'{dir}{self.name}_report.csv',sep=';')
        # save the visualisation
        fig = self.visualise(color_palette)
        fig.savefig(dir+'source_test_hm.png', bbox_inches='tight', dpi=400)


# @TODO
class MultiModelInfoReport(Report):

    def __init__(self) -> None:
        # super().__init__()
        self.table = pd.DataFrame('model','#reactions','#metabolites',
                '#genes','orphans','dead-ends','disconnects','mass unbalanced',
                'charge unbalanced','#reactions with gpr')
        

    def add_single_report(self, report:ModelInfoReport) -> None:
        self.table = pd.concat([self.table,report], ignore_index=True)

    def __add__(self,other):
        self.table = pd.concat(self.table, other.table)

    # @TODO
    def visualise(self):
        pass

    # @TODO
    def save(self):
        pass
