"""Collection of utility functions."""

import libchebipy
import pandas as pd
import requests
from refinegems.analysis_db import BIGG_METABOLITES_URL

__author__ = "Gwendolyn O. Döbel and Carolin Brune"

def add_info_from_ChEBI_BiGG(missing_metabs: pd.DataFrame, charge=True, formula=True, iupac=True) -> pd.DataFrame:
   """Adds information from CHEBI/BiGG to the provided dataframe.
   The following informations can be added:
   - charge
   - formula
   - iupac (name)

   Args:
      - missing_metabs (pd.DataFrame): Table containing metabolites & the respective CHEBI & BiGG IDs
         
   Returns:
      pd.DataFrame: Input table extended with the charges & chemical formulas obtained from CHEBI/BiGG
   """

   # check if a row contains a ChEBI ID, take the first and make sure its in the format: CHEBI:234567
   def get_chebi_id(row: pd.Series) -> str:

      chebi = row.get('ChEBI')
      if pd.isnull(chebi):
         return None
      elif type(chebi) == str: # check for mulitple entries (assuming first one is the best fitting one)
         chebi = chebi.split(',')[0]
         if 'CHEBI' in chebi:
            return chebi
         else:
            return 'CHEBI:' + chebi
      else:
         return 'CHEBI:' + str(chebi)

   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id = get_chebi_id(row)
      bigg_id = str(row.get('bigg_id'))
      charge = None
      if chebi_id:  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         return chebi_entity.get_charge()
      elif bigg_id != 'nan':  # Get charge from BiGG if no ChEBI ID available
         try:
            charge = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['charges'][0]  # Take first charge
         except ValueError:
            pass   
         # If no charge was found, charge=0
         return charge if charge else 0
    
    # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
   def find_formula(row: pd.Series) -> str:
      chebi_id = get_chebi_id(row)
      bigg_id, chem_form = str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id: # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except ValueError:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else 'No formula'
      return chem_formula
   
   # using the ChEBI ID, retrieve the IUPAC name from the ChEBI database.
   def find_iupac(row: pd.Series) -> str:

      chebi_id = get_chebi_id(row)
      if chebi_id:
         chebi_entity = libchebipy.ChebiEntity(chebi_id)
         # only take the IUPAC names
         iupac_names = []
         for name in chebi_entity.get_names():
            print(name)
            if name.get_source() == 'IUPAC':
               iupac_names.append(name.get_name())
               break

         if len(iupac_names) > 0:
            iupac_name = ', '.join(iupac_names)
         else:
            iupac_name = None
      else:
         iupac_name = None

      return iupac_name

   if charge:
      missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   if formula:
      missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
      missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
      missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   if iupac:
      missing_metabs['ChEBI_IUPAC'] = missing_metabs.apply(find_iupac, axis=1)
   
   return missing_metabs