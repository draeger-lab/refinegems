#!/usr/bin/env python
""" Provides functions to compare locus tags (genes) found in BioCyc and in the model

Extracts all BioCyc IDs from the annotations and compares them to a list for your organism from BioCyc.
Reactions with BioCyc IDs not found in the model are expanded to a table containing the BioCyc ID,
the locus tag, the protein ID from NCBI, the EC number, the BiGG ID and the KEGG ID. This section needs 
three TXT files from BioCyc as well as the protein FASTA file from the Genbank entry of the organism. 

These three TXT files can be obtained through creating SmartTables in BioCyc and exporting these as 
Spreadsheets with the parameter FrameID (one & two) or Common name (three). 
SmartTable one: Should contain 'Accession-2' and 'Reactions of gene', 
SmartTable two: Should contain 'Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number', 'Reaction-Direction' and 'Spontaneous?'
SmartTable three: Should contain 'Compound', 'Chemical Formula' and 'InChI-Key'.
"""
# Get all  possible genes by filtering .gff according to 'bio_type=protein_coding' & 'product=hypothetical protein'
# Compare the list of genes with the ones already in the model & add all missing genes
# Before adding to model check if for all genes that are missing for IMITSC147 identifiers exist
# -> Create tables mapping locus tag to old ID, locus tag to new ID & merge 
# -> Specify user input locus_tag start from NCBI PGAP
from libsbml import Model as libModel
import math
import ast
import numpy as np
import pandas as pd
import libchebipy
import requests
from refinegems.entities import get_model_genes, get_model_reacs_or_metabs, compare_gene_lists
from refinegems.analysis_db import get_bigg2other_db, compare_bigg_model, add_stoichiometric_values_to_reacs, BIGG_METABOLITES_URL
from refinegems.io import parse_fasta_headers
import os

__author__ = "Gwendolyn O. Döbel and Dr. Reihaneh Mostolizadeh"


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
def get_biocyc_genes2reactions(inpath: str) -> pd.DataFrame:
   """Parses TSV file from BioCyc to retrieve 'Accession-2' & the corresponding 'Reactions of gene'
   
   Args:
      - inpath (str): Path to file from BioCyc containing the 'Accession-2' to 'Reactions of gene' mapping
      
   Returns:
      pd.DataFrame: Table containing only rows where a 'Reaction of gene' exists
   """
   
   biocyc_genes = pd.read_table(inpath, usecols=['Accession-2', 'Reactions of gene'], dtype=str)
   biocyc_genes.rename(columns={'Accession-2': 'locus_tag', 'Reactions of gene': 'Reaction'}, inplace=True)
   biocyc_genes.replace('', np.nan, inplace=True)
   biocyc_genes.dropna(inplace=True)
   return biocyc_genes


def get_missing_genes2reactions(model_libsbml: libModel, inpath: str) -> pd.DataFrame:
   """Retrieves the missing genes and reactions from the BioCyc table according to the 'Accession-2' identifiers

   Args:
      - model_libsbml (libModel):   Model read in with libSBML
      - inpath (str):               Path to file from BioCyc containing the Accession-2 to Reactions of gene mapping
      
   Returns: 
      pd.DataFrame: Table containing only 'Accession-2' & 'Reactions' for the missing genes
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
   
   Args:
      - inpath (str):   Path to file from BioCyc containing the following columns:
                        'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 
                        'Spontaneous?'
      
   Returns:
      pd.DataFrame: Table containing all biocyc reactions from provided file
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


def extract_metabolites_from_reactions(missing_reactions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
   """Extracts a set of all reactants & products from the missing reactions
   
   Args:
      - missing_reactions (pd.DataFrame): Table containing all missing reactions found through the 
                                       missing genes
                                          
   Returns:
      tuple: Two tables (1) & (2) 
         (1) pd.DataFrame: Table with the column Compound containing all compounds required for the missing BioCyc reactions
         (2) pd.DataFrame: Table with the column bigg_id containing all compounds required for the missing BiGG reactions
   """
   # Get all BioCyc metabolites necessary for the BioCyc reactions
   biocyc_reactants = [r for row in missing_reactions['Reactants'] for r in row]
   biocyc_products = [p for row in missing_reactions['Products'] for p in row]
   biocyc_metabolites = list(set([*biocyc_reactants, *biocyc_products]))
   
   # Get all BiGG metabolites necessary for the BiGG reactions
   bigg_reactions = [ast.literal_eval(str(reac_dict)) for reac_dict in missing_reactions['bigg_reaction']]
   bigg_reactants = [r for reac_dict in bigg_reactions for r in reac_dict.get('reactants')]
   bigg_products = [p for reac_dict in bigg_reactions for p in reac_dict.get('products')]
   bigg_metabolites = list(set([*bigg_reactants, *bigg_products]))
   
   return (pd.DataFrame(biocyc_metabolites, columns=['Compound']), pd.DataFrame(bigg_metabolites, columns=['bigg_id']))


def get_missing_reactions(
   model_libsbml: libModel, genes2reaction: pd.DataFrame, inpath: str
   ) -> tuple[tuple[pd.DataFrame, pd.DataFrame], pd.DataFrame]:
   """Subsets the BioCyc table with the following columns: 
   'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
      to obtain the missing reactions with all the corresponding data 
      & Adds the according BiGG Reaction identifiers

   Args:
      - model_libsbml (libModel):      Model read in with libSBML
      - genes2reaction (pd.DataFrame): Table containing only 'Accession-2' & 'Reactions' for the 
                                       missing genes
      - inpath (str):                  Path to file from BioCyc containing the following columns:
                                       'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 
                                       'Reaction-Direction' 'Spontaneous?'
      
   Returns:
      tuple: Tuple (1) & table (2)
         (1) tuple: Two tables (1) & (2)
            (1) pd.DataFrame: Table containing only the metabolites corresponding to the missing BioCyc reactions
            (2) pd.DataFrame: Table containing only the metabolites corresponding to the missing BiGG reactions
         (2) pd.DataFrame: Table containing the missing reactions with the corresponding data
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
   
   # Add reactants & products dictionary with stoichiometric values to the reactions table
   missing_reactions = add_stoichiometric_values_to_reacs(missing_reactions)
   
   # Get all metabolites for the missing reactions
   biocyc_metabs_from_reacs, bigg_metabs_from_reacs = extract_metabolites_from_reactions(missing_reactions)
   return (biocyc_metabs_from_reacs, bigg_metabs_from_reacs), missing_reactions 


def get_biocyc_metabolites(inpath: str) -> pd.DataFrame:
   """Parses TSV file from BioCyc to retrieve 'Compound (Object ID)' 'Chemical Formula' 'InChI-Key' 'ChEBI'
   
   Args:
      - inpath (str):   Path to file from BioCyc containing the following columns: 
                        'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
      
   Returns:
      pd.DataFrame: Table containing all biocyc metabolites from provided file
   """
   
   biocyc_metabs = pd.read_table(inpath, usecols=['Object ID', 'Chemical Formula', 'InChI-Key', 'ChEBI'], dtype=str)
   biocyc_metabs.rename(columns={'Object ID': 'Compound'}, inplace=True)
   biocyc_metabs.replace('', np.nan, inplace=True)
   biocyc_metabs.dropna(inplace=True)
   return biocyc_metabs

   
def get_missing_metabolites(
   model_libsbml: libModel, metabs_from_reacs: tuple[pd.DataFrame, pd.DataFrame], inpath: str
   ) -> pd.DataFrame: 
   """Subsets the BioCyc table with the following columns: 'Compound' 'Chemical Formula' 'InChI-Key' 'ChEBI'
      to obtain the missing metabolites with all the corresponding data
      & Adds the according BiGG Compound identifiers
      
   Args:
      - model_libsml (libModel):          Model read in with libSBML 
      - metabs_from_reacs (tuple): Two tables containing only the metabolites corresponding to the 
                                          missing reactions for either BioCyc (1) or BiGG (2)
      - inpath (str):                     Path to file from BioCyc containing the following columns:
                                          'Compound' 'Chemical Formula' 'InChI-Key'
      
   Returns:
      tuple: Two tables (1) & (2)
         (1): Table containing the metabolites corresponding to the missing reactions without BiGG IDs
         (2): Table containing the metabolites corresponding to the missing reactions with BiGG IDs
   """
   model_metabs = get_model_reacs_or_metabs(model_libsbml, True)
   biocyc_metabs = get_biocyc_metabolites(inpath)
   biocyc_metabs.rename(columns={'Compound': 'BioCyc'}, inplace=True)
   biocyc_metabs_from_reacs, bigg_metabs_from_reacs = metabs_from_reacs
   
   # Get amount of missing BioCyc metabolites for all missing reactions
   statistics_df.loc['Metabolite', 'Total'] = len(biocyc_metabs_from_reacs['Compound'].unique().tolist())
   
   # Get BiGG BioCyc
   bigg2biocyc_metabs = get_bigg2other_db('BioCyc', True)
   
   # Subset biocyc_metabs with BiGG BioCyc -> To get only metabolites with BiGG IDs
   missing_metabolites = bigg2biocyc_metabs.merge(biocyc_metabs, on='BioCyc') # missing_metabolites
   
   # Filter for all required BiGG metabolites for the missing BiGG reactions
   missing_metabolites = bigg_metabs_from_reacs.merge(missing_metabolites, how='left', on='bigg_id')
   
   # Get amount of missing BioCyc metabolites that have a BiGG ID
   missing_biocyc_metabs = missing_metabolites.dropna(subset=['BioCyc'])
   statistics_df.loc['Metabolite', 'Have BiGG ID'] = len(missing_biocyc_metabs['BioCyc'].unique().tolist())
   
   # Subset missing_metabolites with model_metabs
   missing_metabolites = compare_bigg_model(missing_metabolites, model_metabs, True)
   
   # Get amount of missing metabolites that can & should be added to the model
   statistics_df.loc['Metabolite', 'Can be added'] = len(missing_metabolites['bigg_id'].unique().tolist())
   return missing_metabolites


def get_missing_genes(missing_reactions: pd.DataFrame, fasta: str) -> tuple[pd.DataFrame, pd.DataFrame]:
   """Retrieves all missing genes that belong to the obtained missing reactions
   
   Args: 
      - missing_reactions (pd.DataFrame): Table containing all obtained missing reactions
      - fasta (str):                   Path to a FASTA file where the headers contain the information protein_id and locus_tag
         
   Returns:
      tuple: Two tables (1) & (2)
         (1) pd.DataFrame: Table with the columns locus_tag, Protein_id & Model_id 
                           (The model_id is similar to how CarveMe generates the GeneProduct ID.)
         (2) pd.DataFrame: The input pandas dataframe for the reactions where column 'locus_tag' is exchanged by 'gene_product'
   """
   # Get locus tags from the missing reactions
   locus_tags = list(set([lt for row in missing_reactions['locus_tag'] for lt in row]))
   locus_tags_df = pd.DataFrame(pd.Series(locus_tags), columns=['locus_tag'])
   
   # Get protein and GeneProduct ID for the model from FASTA file
   ids_df = parse_fasta_headers(fasta, id_for_model=True)
   
   # Get the complete dataframe with the protein & model id
   missing_genes = locus_tags_df.merge(ids_df, on='locus_tag')
   statistics_df.loc['Protein', 'Can be added'] = len(missing_genes['locus_tag'].unique().tolist())
   
   # Replace the locus tags in the reaction dataframe with the gene model ID
   def transform_lt_into_gp_model_id(locus_tag_list: list[str]) -> list[str]:
      return [missing_genes.loc[lt, 'model_id'] for lt in locus_tag_list]
   
   missing_genes.set_index('locus_tag', inplace=True)
   missing_reactions['gene_product'] = missing_reactions['locus_tag'].map(transform_lt_into_gp_model_id)
   missing_genes.reset_index(inplace=True)
   missing_reactions.drop('locus_tag', axis=1, inplace=True)
   
   return missing_genes, missing_reactions


def add_charges_chemical_formulae_to_metabs(missing_metabs: pd.DataFrame) -> pd.DataFrame:
   """Adds charges & chemical formulae from CHEBI/BiGG to the provided dataframe

   Args:
      - missing_metabs (pd.DataFrame): Table containing metabolites & the respective CHEBI & BiGG IDs
         
   Returns:
      pd.DataFrame: Input table extended with the charges & chemical formulas obtained from CHEBI/BiGG
   """
   
   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id = str(row.get('bigg_id'))
      charge = None
      if chebi_id:  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
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
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id, chem_form = str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id: # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except (ValueError, IndexError) as e:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else '*'
      return chem_formula
   
   missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
   missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
   missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   
   return missing_metabs

# Inspired by Dr. Reihaneh Mostolizadeh's function to add BioCyc reactions to a model
def replace_reaction_direction_with_fluxes(missing_reacs: pd.DataFrame) -> pd.DataFrame:
   """Extracts the flux lower & upper bounds for each reaction through the entries in column 'Reaction-Direction'
   
   Args:
      - missing_reacs (pd.DataFrame): Table containing reactions & the respective Reaction-Directions
         
   Returns:
      pd.DataFrame: Input table extended with the fluxes lower & upper bounds obtained from 
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
   model_libsbml: libModel, biocyc_file_paths: list[str]
   ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
   """Main function to retrieve the tables for missing genes, metabolites and reactions from BioCyc
   
   Args:
      - model_libsbml (libModel):   libSBML Model object
      - biocyc_file_paths (list):   List of the files required for the BioCyc analysis
         
   Returns: 
      tuple: Five tables (1) - (5)
         (1) pd.DataFrame: Table containing the statistics of the BioCyc gapfill analysis
         (2) pd.DataFrame: Table containing the missing genes that belong to the missing reactions
         (3) pd.DataFrame: Table containing the missing metabolites with BiGG IDs belonging to the missing reactions
         (4) pd.DataFrame: Table containing the missing metabolites without BiGG IDs belonging to the missing reactions
         (5) pd.DataFrame: Table containing the missing reactions
   """
   # check paths
   for path in biocyc_file_paths:
      if not os.path.exists(path):
         print(f"{path} does not exists, check input!")
         break

   # Extract missing reactions from all missing genes
   genes2reactions = get_missing_genes2reactions(model_libsbml, biocyc_file_paths[0])
   metabs_from_reacs, missing_reactions_df = get_missing_reactions(model_libsbml, genes2reactions, biocyc_file_paths[1])
   missing_reactions_df = replace_reaction_direction_with_fluxes(missing_reactions_df)
   
   # Extract missing metabolites that belong to the missing reactions
   missing_metabolites_df= get_missing_metabolites(model_libsbml, metabs_from_reacs, biocyc_file_paths[2]) 
   missing_metabolites_df = add_charges_chemical_formulae_to_metabs(missing_metabolites_df)
   
   # Extract missing genes that belong to the missing reactions
   missing_genes_df, missing_reactions_df = get_missing_genes(missing_reactions_df, biocyc_file_paths[3])
   
   # Remove index from statistics_df
   statistics_df.reset_index(inplace=True)
   
   return (statistics_df, missing_genes_df, missing_metabolites_df, missing_reactions_df) 
