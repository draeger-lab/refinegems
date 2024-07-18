"""_summary_
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

############################################################################
# requirements
############################################################################

from abc import ABC, abstractmethod

############################################################################
# variables
############################################################################

############################################################################
# classes
############################################################################

class GapFiller(ABC):

   def __init__(self) -> None:
      self.full_gene_list = None
      self._statistics = dict()
    
   @abstractmethod
   def load_gene_list(self):
      pass

   def get_missing_genes(self, model):

      geneps_in_model = [_.id for _ in model.genes]

      return [_ for _ in self.full_gene_list if _ not in geneps_in_model]

   def fill_model(self, model):
      pass

   def calculate_stats(self):
      pass

   def report(self):
      pass


class KEGGapFiller(GapFiller):
   pass

class BioCycGapFiller(GapFiller):
   pass

class GeneGapFiller(GapFiller):
   pass

############################################################################
# functions
############################################################################