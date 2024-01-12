"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune, Famke Baeuerle, Gwendolyn O. Döbel'

################################################################################
# requirements
################################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings

from importlib.resources import files
from pathlib import Path
from typing import Literal

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
    pass


class SingleGrowthSimulationReport(Report):
    
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

    def __init__(self, reports = []):

        self.reports = reports
        self.models = set([_.model_name for _ in reports])
        self.media = set([_.medium_name for _ in reports])

    def __str__(self):
        
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
                cmap = matplotlib.cm.get_cmap(color_palette).copy()
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.cm.get_cmap('YlGn').copy()

            # set up the figure
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])

            # construct the plot
            max_ydata = max(ydata)
            if max_ydata <= 0:
                warnings.warn('Model is not able to grow on every medium. Returning empty figure.')
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
                cmap = matplotlib.cm.get_cmap(color_palette).copy()
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.cm.get_cmap('YlGn').copy()
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
            cmap = matplotlib.cm.get_cmap(color_palette).copy()
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.cm.get_cmap('YlGn').copy()
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
            cmap = matplotlib.cm.get_cmap(color_palette).copy()
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.cm.get_cmap('YlGn').copy()

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
        
        
    def save(self, dir:str, color_palette:str='YnGr'):
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