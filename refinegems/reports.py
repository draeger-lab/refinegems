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
from typing import Literal

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

        