"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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

    def add_sim_results(self, new_rep: SingleGrowthSimulationReport):
        """Add a new single growth report to the reports list

        Args:
            new_rep (SingleGrowthSimulationReport): The new simulation report.
        """

        self.reports.append(new_rep)

    def to_table(self):

        l = []
        for report in self.reports:
            l.append(report.to_dict())
        return pd.DataFrame(l)

    # @REWRITE to make it work with the new standard
    # --> only heatmap if n-m
    def plot_heatmap_dt(self):
        """Creates heatmap of simulated doubling times with additives
        
        Args:
            - growth (pd.DataFrame): Containing growth data from simulate_all
            
        Returns:
            plot: Seaborn Heatmap
        """
        growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
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
        cmap=matplotlib.cm.get_cmap('YlGn').copy()
        cmap.set_under('black')
        cmap.set_over('white')
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
    

    def plot_heatmap_native(growth: pd.DataFrame):
        """Creates a plot were if growth without additives is possible is marked from yellow to green otherwise black

        Args:
            - growth (pd.DataFrame): Containing growth data from simulate_all
            
        Returns:
            plot: Seaborn Heatmap
        """
        def get_native_growth(row):
            if row['complete'] == True:
                return row['doubling_time [min]']
            else:
                return 0
        
        growth['native_growth'] = growth.apply(get_native_growth, axis=1)
        growth = growth[['medium', 'model', 'native_growth']]
        growth=growth.set_index(['medium', 'model']).sort_index().T.stack()
        growth.columns.name=None
        growth.index.names = (None,None)
        growth.index.name=None
        growth.index = growth.index.get_level_values(1)
        growth[growth > 500] = 0
        growth[growth < 0] = 0
        growth.replace([np.inf, -np.inf], 0, inplace=True)
        over_growth = growth.max().max() + 6
        growth.replace(np.nan, over_growth, inplace=True)
        annot = growth.copy()
        annot = annot.round().astype(int)
        annot[annot == np.nan] = 'No data'
        annot[annot < 1e-5] = ''
        annot.replace(over_growth.round().astype(int), 'No data', inplace=True)
        under_growth = growth.min().min() - 5
        vmin= under_growth if under_growth > 1e-5 else 1e-5 #Use same threshhold as in find_missing_essential in growth
        vmax= over_growth - 1
        cmap=matplotlib.cm.get_cmap('YlGn').copy()
        cmap.set_under('black')
        cmap.set_over('white')
        fig, ax = plt.subplots(figsize=(10,8))
        sns.heatmap(growth.T,
                    annot=annot.T, 
                    annot_kws={"fontsize":15},
                    vmin=vmin,
                    vmax=vmax,
                    cmap=cmap,
                    linewidth=.5, 
                    ax=ax,
                    cbar_kws={'orientation':'vertical', 'label':'doubling time [min]', 'extend': 'min', 'extendrect':True},
                    fmt=''
                    )
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)
        rotation = 40 if len(growth.index) > 3 else 0
        plt.tick_params(rotation=0, bottom=False, top=False, left=False, right=False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha="right")
        return fig