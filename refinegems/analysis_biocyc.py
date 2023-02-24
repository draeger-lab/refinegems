#!/usr/bin/env python
""" Provides functions to compare locus tags (genes) found in BioCyc and in the model

Extracts all BioCyc IDs from the annotations and compares them to a list for your organism from BioCyc.
Reactions with BioCyc IDs not found in the model are expanded to a table containing the BioCyc ID,
the locus tag, the protein ID from NCBI, the EC number, the BiGG ID and the KEGG ID.
This section needs the Genbank GFF file of your organism, the TXT file from BiGG containing all reactions as well as 
the TXT file from BiGG containing all metabolites and three TXT files from BioCyc. These three TXT files can be 
obtained through creating SmartTables in BioCyc and exporting these as Spreadsheets with the parameter FrameID.
SmartTable one should contain 'Accession-2' and 'Reactions of gene', SmartTable two should contain 'Reaction', 
'Reactants of reaction', 'Products of reaction', 'EC-Number', 'Reaction-Direction' and 'Spontaneous?' while SmartTable 
three should contain 'Compound', 'Chemical Formula' and 'InChI-Key'.

Due to the KEGG REST API this is relatively slow (model of size 1500 reactions - 20 min).
"""
# Get all  possible genes by filtering .gff according to 'bio_type=protein_coding' & 'product=hypothetical protein'
# Compare the list of genes with the ones already in the model & add all missing genes
# Before adding to model check if for all genes that are missing for IMITSC147 identifiers exist
# -> Create dataframes mapping locus tag to old ID, locus tag to new ID & merge 
# -> Specify user input locus_tag start from NCBI PGAP
from libsbml import *
import numpy as np
import pandas as pd
import libchebipy
import requests
from refinegems.entities import get_model_genes, get_model_reacs_or_metabs, compare_gene_lists
from refinegems.analysis_db import get_bigg2other_db, compare_bigg_model, add_stoichiometric_values_to_reacs
from refinegems.io import parse_fasta_headers

__author__ = "Gwendolyn O. Gusak"


# Global variable for statistics
statistics_dict = {
   'Missing entity': ['Protein', 'Metabolite', 'Reaction'], 
   'Total': [np.NaN, np.NaN, np.NaN], 
   'Have BiGG ID': [np.NaN, np.NaN, np.NaN], 
   'Can be added': [np.NaN, np.NaN, np.NaN],
   'Notes': [
      'Amount derived from locus tag comparison/Amount remaining for the reactions',
      'Only metabolites that are required for the missing reactions',
      'Only reactions that belong to the missing genes/proteins'
      ]
   }
statistics_df = pd.DataFrame(statistics_dict).set_index('Missing entity')


# Locus tags in GenBank GFF file == BioCyc Accession-2 == Old locus tags in RefSeq GFF file
# Locus tags in RefSeq GFF file == BioCyc Accession-1
# Label in model == Locus tag from GenBank GFF file == BioCyc Accession-2
def get_biocyc_genes2reactions(inpath: str):
   """Parses TSV file from BioCyc to retrieve 'Accession-2' & the corresponding 'Reactions of gene'
   
      Params:
         - inpath (str): Path to file from BioCyc containing the 'Accession-2' to 'Reactions of gene' mapping
      
      Returns:
         -> Pandas dataframe containing only rows where a 'Reaction of gene' exists
   """
   
   biocyc_genes = pd.read_table(inpath, usecols=['Accession-2', 'Reactions of gene'], dtype=str)
   biocyc_genes.rename(columns={'Accession-2': 'locus_tag', 'Reactions of gene': 'Reaction'}, inplace=True)
   biocyc_genes.replace('', np.nan, inplace=True)
   biocyc_genes.dropna(inplace=True)
   return biocyc_genes


def get_missing_genes2reactions(model_libsbml: Model, inpath:str) -> tuple[list[str], pd.DataFrame]:
   """Retrieves the missing genes and reactions from the BioCyc table according to the 'Accession-2' identifiers

      Params:
         - model_libsbml (Model):   Model read in with libSBML
         - inpath (str):            Path to file from BioCyc containing the Accession-2 to Reactions of gene mapping
      
      Returns:
         -> A pandas dataframe containing only 'Accession-2' & 'Reactions' for the missing genes
   """
   
   gps_in_model = get_model_genes(model_libsbml)
   biocyc_genes = get_biocyc_genes2reactions(inpath)
   missing_biocyc_genes = compare_gene_lists(gps_in_model, biocyc_genes, False)
   
   missing_biocyc_reactions = pd.DataFrame(
      missing_biocyc_genes['Reaction'].str.split('//').tolist(), index=missing_biocyc_genes['locus_tag']
      ).stack()
   missing_biocyc_reactions = missing_biocyc_reactions.reset_index([0, 'locus_tag'])
   missing_biocyc_reactions.columns = ['locus_tag', 'Reaction']
   missing_biocyc_reactions['Reaction'] = missing_biocyc_reactions['Reaction'].str.strip()
   
   # Get amount of missing genes from BioCyc for statistics
   biocyc_genes = missing_biocyc_genes['locus_tag'].unique().tolist()
   statistics_df.loc['Protein', 'Total'] = len(biocyc_genes)
   
   return missing_biocyc_reactions


def get_biocyc_reactions(inpath: str) -> pd.DataFrame:
   """Parses TSV file from BioCyc to retrieve 'Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number',
      'Reaction-Direction' & 'Spontaneous?'
   
      Params:
         - inpath (str):   Path to file from BioCyc containing the following columns:
                           'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 
                           'Spontaneous?'
      
      Returns:
         -> A pandas dataframe containing all biocyc reactions from provided file
   """
   
   biocyc_reacs = pd.read_table(inpath, usecols=
                                ['Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number', 
                                 'Reaction-Direction', 'Spontaneous?'],
                                dtype=str
                                )
   biocyc_reacs.rename(columns=
                       {'Reactants of reaction': 'Reactants', 'Products of reaction': 'Products', 'EC-Number': 'EC'},
                       inplace=True
                       )
   biocyc_reacs.replace('', np.nan, inplace=True)
   biocyc_reacs['Spontaneous?'] = biocyc_reacs['Spontaneous?'].fillna('F')
   biocyc_reacs.dropna(subset=['Reaction', 'Reactants', 'Products', 'EC'], inplace=True)
   return biocyc_reacs


def extract_metabolites_from_reactions(missing_reactions: pd.DataFrame):
   """Extracts a set of all reactants & products from the missing reactions
   
      Params:
         - missing_reactions (DataFrame): A pandas dataframe containing all missing reactions found through the 
                                          missing genes
                                          
      Returns:
         -> A pandas dataframe with the column Compound containing all compounds required for the missing reactions
   """
   reactants = [r for row in missing_reactions['Reactants'] for r in row]
   products = [p for row in missing_reactions['Products'] for p in row]
   metabolites = list(set([*reactants, *products]))
   return pd.DataFrame(metabolites, columns=['Compound'])


def get_missing_reactions(
   model_libsbml: Model, genes2reaction: pd.DataFrame, inpath: str
   ) -> tuple[pd.DataFrame, pd.DataFrame]:
   """Subsets the BioCyc table with the following columns: 
   'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
      to obtain the missing reactions with all the corresponding data 
      & Adds the according BiGG Reaction identifiers

      Params:
         - model_libsbml (Model):         Model read in with libSBML
         - genes2reaction (pd.DataFrame): A pandas dataframe containing only 'Accession-2' & 'Reactions' for the 
                                          missing genes
         - inpath (str):                  Path to file from BioCyc containing the following columns:
                                          'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 
                                          'Reaction-Direction' 'Spontaneous?'
      
      Returns:
         -> Two pandas dataframes (1) & (2)
            (1): A pandas dataframe containing only the metabolites corresponding to the missing reactions
            (2): A pandas dataframe containing the missing reactions with the corresponding data
   """
   model_reacs = get_model_reacs_or_metabs(model_libsbml)
   biocyc_reacs = get_biocyc_reactions(inpath)
   
   # Get missing reactions from missing genes
   missing_reactions = genes2reaction.merge(biocyc_reacs, on='Reaction')
   
   # Turn entries with '//' into lists
   missing_reactions['Reactants'] = missing_reactions['Reactants'].str.split('\s*//\s*')
   missing_reactions['Products'] = missing_reactions['Products'].str.split('\s*//\s*')
   missing_reactions['EC'] = missing_reactions['EC'].str.split('\s*//\s*')
   
   # Turn locus_tag column into lists of locus tags per reaction
   locus_tags_as_list = missing_reactions.groupby('Reaction')['locus_tag'].apply(list).reset_index(name='locus_tag')
   missing_reactions.drop('locus_tag', axis=1, inplace=True)
   missing_reactions = locus_tags_as_list.merge(missing_reactions, on='Reaction')
   statistics_df.loc['Reaction', 'Total'] = len(missing_reactions['Reaction'].unique().tolist())
   
   # Get BiGG BioCyc
   bigg2biocyc_reacs = get_bigg2other_db('BioCyc')
   
   # Subset missing_reactions with BiGG BioCyc
   missing_reactions.rename(columns={'Reaction': 'BioCyc'}, inplace=True)
   missing_reactions = bigg2biocyc_reacs.merge(missing_reactions, on='BioCyc')
   
   # Get amount of missing reactions that have a BiGG ID
   statistics_df.loc['Reaction', 'Have BiGG ID'] = len(missing_reactions['BioCyc'].unique().tolist())
   
   # Subset missing_reactions with model_reacs
   missing_reactions = compare_bigg_model(missing_reactions, model_reacs)
   
   # Get amount of missing reactions that are not in the model
   statistics_df.loc['Reaction', 'Can be added'] = len(missing_reactions['bigg_id'].unique().tolist())
   
   # Get all metabolites for the missing reactions
   metabs_from_reaction = extract_metabolites_from_reactions(missing_reactions)
   return metabs_from_reaction, missing_reactions  


def get_biocyc_metabolites(inpath: str) -> pd.DataFrame:
   """Parses TSV file from BioCyc to retrieve 'Compound (Object ID)' 'Chemical Formula' 'InChI-Key' 'ChEBI'
   
      Params:
         - inpath (str):   Path to file from BioCyc containing the following columns: 
                           'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
      
      Returns:
         -> A pandas dataframe containing all biocyc metabolites from provided file
   """
   
   biocyc_metabs = pd.read_table(inpath, usecols=['Object ID', 'Chemical Formula', 'InChI-Key', 'ChEBI'], dtype=str)
   biocyc_metabs.rename(columns={'Object ID': 'Compound'}, inplace=True)
   biocyc_metabs.replace('', np.nan, inplace=True)
   biocyc_metabs.dropna(inplace=True)
   return biocyc_metabs


def get_missing_metabolites_wo_BiGG(
   missing_metabs_overall: pd.DataFrame, missing_metabs_BiGG: pd.DataFrame
   ) -> pd.DataFrame:
   """Retrieves all missing metabolites that have no BiGG ID mapping
   
      Params:
         - missing_metabs_overall (DataFrame):  A pandas dataframe containing all missing metabolites retrieved from
                                                all missing reactions
         - missing_metabs_BiGG (DataFrame):     A pandas dataframe containing all missing metabolites with BiGG ID mappings
         
      Returns:
         -> A pandas dataframe containing all missing metabolites without BiGG ID mappings
   """
   all_missing_metabs = missing_metabs_overall.set_index('BioCyc')
   metabs_with_BiGG = missing_metabs_BiGG.set_index('BioCyc')
   
   remaining_missing_metabs = all_missing_metabs[~all_missing_metabs.index.isin(metabs_with_BiGG.index)]
    
   return remaining_missing_metabs.reset_index()
   

def get_missing_metabolites(
   model_libsbml: Model, metabs_from_reacs: pd.DataFrame, inpath: str
   ) -> tuple[pd.DataFrame, pd.DataFrame]:
   """Subsets the BioCyc table with the following columns: 'Compound' 'Chemical Formula' 'InChI-Key' 'ChEBI'
      to obtain the missing metabolites with all the corresponding data
      & Adds the according BiGG Compound identifiers
      
      Params:
         - model_libsml (Model):             Model read in with libSBML 
         - metabs_from_reacs (pd.DataFrame): A pandas dataframe containing only the metabolites corresponding to the 
                                             missing reactions
         - inpath (str):                     Path to file from BioCyc containing the following columns:
                                             'Compound' 'Chemical Formula' 'InChI-Key'
      
      Returns:
         -> Two pandas dataframes (1) & (2)
            (1): A pandas dataframe containing the metabolites corresponding to the missing reactions without BiGG IDs
            (2): A pandas dataframe containing the metabolites corresponding to the missing reactions with BiGG IDs
   """
   model_metabs = get_model_reacs_or_metabs(model_libsbml, True)
   biocyc_metabs = get_biocyc_metabolites(inpath)
   
   # Get missing metabolites for missing reactions
   missing_metabolites = metabs_from_reacs.merge(biocyc_metabs, how='left', on='Compound')
   missing_metabolites.rename(columns={'Compound': 'BioCyc'}, inplace=True)
   missing_metabs_overall = missing_metabolites
   
   # Get amount of missing metabolites for all missing reactions
   statistics_df.loc['Metabolite', 'Total'] = len(missing_metabolites['BioCyc'].unique().tolist())
   
   # Get BiGG BioCyc
   bigg2biocyc_metabs = get_bigg2other_db('BioCyc', True)
   
   # Subset missing_metabolites with BiGG BioCyc -> To get only metabolites with BiGG IDs
   missing_metabolites = bigg2biocyc_metabs.merge(missing_metabolites, on='BioCyc')
   
   # Retrieve missing metabolites that have no BiGG IDs but belong to the missing reactions
   missing_metabs_wo_BiGG = get_missing_metabolites_wo_BiGG(missing_metabs_overall, missing_metabolites)
   
   # Get amount of missing metabolites that have a BiGG ID
   statistics_df.loc['Metabolite', 'Have BiGG ID'] = len(missing_metabolites['BioCyc'].unique().tolist())
   
   # Subset missing_metabolites with model_metabs
   missing_metabolites = compare_bigg_model(missing_metabolites, model_metabs)
   
   # Get amount of missing metabolites that are not in the model
   statistics_df.loc['Metabolite', 'Can be added'] = len(missing_metabolites['bigg_id'].unique().tolist())
   return missing_metabolites, missing_metabs_wo_BiGG


def get_missing_genes(missing_reactions: pd.DataFrame, fasta: str) -> pd.DataFrame:
   """Retrieves all missing genes that belong to the obtained missing reactions
   
      Params: 
         - missing_reactions (DataFrame): A pandas dataframe containing all obtained missing reactions
         - fasta (str):                   Path to a FASTA file where the headers contain the information protein_id and locus_tag
         
      Returns:
         -> A pandas dataframe with the columns locus_tag, Protein_id & Model_id 
         (The model_id is similar to how CarveMe generates the GeneProduct ID.)
   """
   # Get locus tags from the missing reactions
   locus_tags = list(set([lt for row in missing_reactions['locus_tag'] for lt in row]))
   locus_tags_df = pd.DataFrame(pd.Series(locus_tags), columns=['locus_tag'])
   
   # Get protein and GeneProduct ID for the model from FASTA file
   ids_df = parse_fasta_headers(fasta, id_for_model=True)
   
   # Get the complete dataframe with the protein & model id
   missing_genes = locus_tags_df.merge(ids_df, on='locus_tag')
   statistics_df.loc['Protein', 'Can be added'] = len(missing_genes['locus_tag'].unique().tolist())
   
   return missing_genes


def add_charges_chemical_formulae_to_metabs(missing_metabs: pd.DataFrame) -> pd.DataFrame:
   """Adds charges & chemical formulae from CHEBI to the provided dataframe

      Params:
         - missing_metabs (DataFrame): A pandas dataframe containing metabolites & the respective CHEBI IDs
         
      Returns:
         -> The input pandas dataframe extended with the charges & chemical formulas obtained from CHEBI
   """
   metab_bigg_url = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/'
   
   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id, bigg_id = str(row.get('ChEBI')), str(row.get('bigg_id'))
      charge = None
      if chebi_id != 'nan':  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         return chebi_entity.get_charge()
      elif bigg_id != 'nan':  # Get charge from BiGG if no ChEBI ID available
         try:
            charge = requests.get(metab_bigg_url + bigg_id[:-2]).json()['charges'][0]  # Take first charge
         except ValueError:
            pass   
         # If no charge was found, charge=0
         return charge if charge else 0
   
   # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
   def find_formula(row: pd.Series) -> str:
      chebi_id, bigg_id, chem_form = str(row.get('ChEBI')), str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id != 'nan': # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(metab_bigg_url + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except ValueError:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else 'No formula'
      return chem_formula
   
   missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
   missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
   missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   
   return missing_metabs

# Inspired by Dr. Reihaneh Mostolizadeh's function to add BioCyc reactions to a model
def replace_reaction_direction_with_fluxes(missing_reacs: pd.DataFrame) -> pd.DataFrame:
   """Extracts the flux lower & upper bounds for each reaction through the entries in column 'Reaction-Direction'
   
      Params:
         - missing_reacs (DataFrame): A pandas dataframe containing reactions & the respective Reaction-Directions
         
      Returns:
         -> The input pandas dataframe extended with the fluxes lower & upper bounds obtained from 
            the Reaction-Directions
   """
    
   def get_fluxes(row: pd.Series) -> dict[str: str]:
      direction = row['Reaction-Direction']
      fluxes = {}
      
      if type(direction) == float:
         # Use default bounds as described in readthedocs from COBRApy
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'RIGHT-TO-LEFT' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_0_bound'
      elif 'LEFT-TO-RIGHT' in direction:
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'REVERSIBLE' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_default_ub'
      
      return str(fluxes)
   
   missing_reacs['fluxes'] = missing_reacs.apply(get_fluxes, axis=1)
   missing_reacs.drop('Reaction-Direction', axis=1, inplace=True)
   
   return missing_reacs


def biocyc_gene_comp(
   model_libsbml: Model, biocyc_file_paths: list[str]
   ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
   """Main function to retrieve the dataframes for missing genes, metabolites and reactions from BioCyc
   
      Params:
         - model_libsbml (Model):      libSBML Model object
         - biocyc_file_paths (list):   A list of the files required for the BioCyc analysis
         
      Returns: 
         -> Five dataframes (1) - (5):
            (1): A pandas dataframe containing the statistics of the BioCyc gapfill analysis
            (2): A pandas dataframe containing the missing genes that belong to the missing reactions
            (3): A pandas dataframe containing the missing metabolites with BiGG IDs belonging to the missing reactions
            (4): A pandas datafrane containing the missing metabolites without BiGG IDs belonging to the missing reactions
            (5): A pandas dataframe containing the missing reactions
   """
   # Extract missing reactions from all missing genes
   genes2reactions = get_missing_genes2reactions(model_libsbml, biocyc_file_paths[0])
   metabs_from_reacs, missing_reactions_df = get_missing_reactions(model_libsbml, genes2reactions, biocyc_file_paths[1])
   missing_reactions_df = add_stoichiometric_values_to_reacs(missing_reactions_df)
   missing_reactions_df = replace_reaction_direction_with_fluxes(missing_reactions_df)
   
   # Extract missing metabolites that belong to the missing reactions
   missing_metabolites_df, missing_metabs_wo_BiGG_df = get_missing_metabolites(model_libsbml, metabs_from_reacs, biocyc_file_paths[2])
   missing_metabolites_df = add_charges_chemical_formulae_to_metabs(missing_metabolites_df)
   # If metabolites should be added due to reactions but no BiGG ID was found
   if len(missing_metabs_wo_BiGG_df.index) > 0:
      missing_metabs_wo_BiGG_df = add_charges_chemical_formulae_to_metabs(missing_metabs_wo_BiGG_df)
   
   # Extract missing genes that belong to the missing reactions
   missing_genes_df = get_missing_genes(missing_reactions_df, biocyc_file_paths[3])
   
   # Remove index from statistics_df
   statistics_df.reset_index(inplace=True)
   
   return (statistics_df, missing_genes_df, missing_metabolites_df, missing_metabs_wo_BiGG_df, missing_reactions_df)
   