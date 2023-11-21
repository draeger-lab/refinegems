"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################


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


class GrowthSimulationReport(Report):

    def __init__(self, reports = []):

        self.reports = reports

    def add_sim_results(self, new_rep: SingleGrowthSimulationReport):
        
        self.reports.append(new_rep)