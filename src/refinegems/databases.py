#!/usr/bin/env python
import io
import re
import sqlite3
import requests
import pandas as pd
from enum import Enum
from sqlite3 import Error
from os import path
from importlib.resources import files

__author__ = 'Gwendolyn O. Döbel'

PATH_TO_DB = files('refinegems.data.database').joinpath('data.db')
VERSION_FILE = files('refinegems.data.database').joinpath('current_bigg_db_version.txt') 
VERSION_URL = 'http://bigg.ucsd.edu/api/v2/database_version'

class ValidationCodes(Enum):
   """Validation codes for the database

      Args:
         - Enum (Enum): Provided as input to get a number mapping for the codes
   """
   COMPLETE = 0,  # All tables are in data.db
   EMPTY = 1,  # data.db is either empty or incorrect
   BIGG = 2,  # Only BiGG tables are in data.db
   SBO_MEDIA = 3,  # Only SBO & Media tables are in data.db (Can only occurr together) 
   BIGG_SBO_MEDIA = 4,  # Only BiGG, SBO and media tables are in data.db 
   MODELSEED_COMPOUNDS = 5,  # Only ModelSEED compounds table is in data.db
   BIGG_MSEED_COMPPOUNDS = 6,  # Only Bigg and ModelSEED compounds tables are in data.db
   SBO_MEDIA_MSEED_COMPOUNDS = 7  # Only SBO, media and ModelSEED compounds tables are in data.db


def is_valid_database(db_cursor: sqlite3.Cursor) -> int:
   """Verifies if database has:
         - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
         - 2 tables with names 'bigg_to_sbo' & 'ec_to_sbo'
         - 2 tables with names 'media' & 'media_composition'
         - 1 table with name 'modelseed_compounds'
   
   Args:
      - db_cursor (sqlite3.Cursor): Cursor from open connection to the database (data.db)

   Returns:
      int: Corresponding to one of the ValidationCodes
   """
   print('Verifying database...')
   
   # Fetches the table names as string tuples from the connected database
   db_cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
   tables = [string[0] for string in db_cursor.fetchall()]
   
   bigg_tables_contained = len([s for s in tables if re.match('^bigg_(?!to)(.*?)', s, re.IGNORECASE)]) == 2
   sbo_tables_contained = len([s for s in tables if re.match('(.*?)_sbo$', s, re.IGNORECASE)]) == 2
   media_tables_contained = len([s for s in tables if re.match('media', s, re.IGNORECASE)]) == 2
   sbo_media_tables_contained = sbo_tables_contained and media_tables_contained  # These can only occur together
   modelseed_cmpd_tbl_contained = len([s for s in tables if s == 'modelseed_compounds']) == 1
   
   bigg_sbo_media_tbls_contained = bigg_tables_contained and sbo_media_tables_contained
   bigg_modelseed_cmpd_tbls_contained = bigg_tables_contained and modelseed_cmpd_tbl_contained
   sbo_media_modelseed_cmpd_tbls_contained = sbo_media_tables_contained and modelseed_cmpd_tbl_contained
   all_tables_contained = bigg_sbo_media_tbls_contained and modelseed_cmpd_tbl_contained
   
   if all_tables_contained: return ValidationCodes.COMPLETE
   elif bigg_modelseed_cmpd_tbls_contained: return ValidationCodes.BIGG_MSEED_COMPOUNDS
   elif sbo_media_modelseed_cmpd_tbls_contained: return ValidationCodes.SBO_MEDIA_MSEED_COMPOUNDS
   elif bigg_sbo_media_tbls_contained: return ValidationCodes.BIGG_SBO_MEDIA
   elif bigg_tables_contained: return ValidationCodes.BIGG
   elif sbo_media_tables_contained: return ValidationCodes.SBO_MEDIA
   elif modelseed_cmpd_tbl_contained: return ValidationCodes.MODELSEED_COMPOUNDS
   else: return ValidationCodes.EMPTY


def create_sbo_media_database(db_cursor: sqlite3.Cursor):
   """Creates the SBO annotation database with 2 tables ('bigg_to_sbo' & 'ec_to_sbo')
      & the media database with 2 tables ('media', 'media_compositions') from file './data/database/sbo_media_db.sql'

   Args:
      - db_cursor (sqlite3.Cursor): Cursor from open connection to the database (data.db)
   """
   print('Adding SBO and media tables...')
   
   with open(files('refinegems.data.database').joinpath('sbo_media_db.sql')) as schema:
      db_cursor.executescript(schema.read())


def update_bigg_db(latest_version: str, db_connection: sqlite3.Connection):
   """Updates the BiGG tables 'bigg_metabolites' & 'bigg_reactions' within a database (data.db)

   Args:
      - latest_version (str): String containing the latest version of the BiGG database
      - db_connection (sqlite3.Connection): Open connection to the database (data.db)
   """
   print('Adding BiGG tables...')
   db_connection.execute('DROP TABLE IF EXISTS bigg_metabolites')
   db_connection.execute('DROP TABLE IF EXISTS bigg_reactions')
   
   # Store currently used version
   with open(VERSION_FILE, 'w') as file:
      file.write(latest_version)
   
   # Create BiGG metabolites table
   BIGG_MODELS_METABS_URL = 'http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt'
   bigg_models_metabs = requests.get(BIGG_MODELS_METABS_URL).text
   bigg_models_metabs_df = pd.read_csv(io.StringIO(bigg_models_metabs), dtype=str, sep='\t')
   bigg_models_metabs_df.to_sql('bigg_metabolites', db_connection, index=False)

   # Create BiGG reactions table
   BIGG_MODELS_REACS_URL = 'http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
   bigg_models_reacs = requests.get(BIGG_MODELS_REACS_URL).text
   bigg_models_reacs_df = pd.read_csv(io.StringIO(bigg_models_reacs), dtype=str, sep='\t')
   bigg_models_reacs_df.to_sql('bigg_reactions', db_connection, index=False)
   

def get_latest_bigg_databases(db_connection: sqlite3.Connection, is_missing: bool=True):
   """Gets the latest BiGG tables for metabolites & reactions if:
         - No version file is locally available
         - The version in the local version file is NOT the latest
         - No BiGG tables currently exist in the database

   Args:
      - db_connection (sqlite3.Connection): Open connection to the database (data.db)
      - is_missing (bool, optional): True if no BiGG tables are in the database. Defaults to True.
   """
   # Check if BiGG database had an update
   LATEST_VERSION = requests.get(VERSION_URL).json()['bigg_models_version']
   
   if not path.exists(VERSION_FILE) or is_missing:
      update_bigg_db(LATEST_VERSION, db_connection)
      
   else:
      with open(VERSION_FILE, 'r') as file:
         version = file.readline().strip()
         
      if version != LATEST_VERSION:
         update_bigg_db(LATEST_VERSION, db_connection)
 
 
def get_modelseed_compounds_database(db_connection: sqlite3.Connection):
   """Retrieves the compounds table from ModelSEED from the respective GitHub repository

   Args:
      - db_connection (sqlite3.Connection): Open connection to the database (data.db)
   """
   print('Adding the ModelSEED compounds table...')
   MODELSEED_COMPOUNDS_URL = 'https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv'
   modelseed_compounds = requests.get(MODELSEED_COMPOUNDS_URL).text
   modelseed_df = pd.read_csv(io.StringIO(modelseed_compounds), sep='\t')
   modelseed_df.to_sql('modelseed_compounds', db_connection, index=False, if_exists='replace')
    
         
def initialise_database():
   """Initialises/updates the database (data.db)

      After initialisation the database contains:
         - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
         - 2 tables with names 'bigg_to_sbo' & 'ec_to_sbo'
         - 2 tables with names 'media' & 'media_composition'
         - 1 table with name 'modelseed_compounds' 
   """
   # Initialise empty connection
   con = None

   print('Initialising database...')
   
   # Try to open connection & get cursor
   try:
      con = sqlite3.connect(PATH_TO_DB)
      cursor = con.cursor()
      
      validity_code = is_valid_database(cursor)
      
      if validity_code == ValidationCodes.BIGG:
         print('Only BiGG tables contained in database.')
         create_sbo_media_database(cursor)
         get_modelseed_compounds_database(con)
      
      elif validity_code == ValidationCodes.SBO_MEDIA:
         print('Only SBO and media tables contained in database.')
         get_latest_bigg_databases(con)
         get_modelseed_compounds_database(con)
         
      elif validity_code == ValidationCodes.MODELSEED_COMPOUNDS:
         print('Only ModelSEED compounds table contained in database.')
         create_sbo_media_database(cursor)
         get_latest_bigg_databases(con)
         
      elif validity_code == ValidationCodes.BIGG_SBO_MEDIA:
         print('Only BiGG, SBO and media tables contained in database.')
         get_modelseed_compounds_database(con)
         
      elif validity_code == ValidationCodes.BIGG_MSEED_COMPPOUNDS:
         print('Only BiGG and ModelSEED compounds tables contained in database.')
         create_sbo_media_database(cursor)
         
      elif validity_code == ValidationCodes.SBO_MEDIA_MSEED_COMPOUNDS:
         print('Only SBO, media and ModelSEED compounds tables contained in database.')
         get_latest_bigg_databases(con)
         
      elif validity_code == ValidationCodes.EMPTY:
         print('Incorrect or empty database. Initialise database with required tables...')
         create_sbo_media_database(cursor)
         get_latest_bigg_databases(con)
         get_modelseed_compounds_database(con)
         
      elif validity_code == ValidationCodes.COMPLETE:
         print('Verifying if BiGG tables are up-to-date...')
         get_latest_bigg_databases(con, False)
      
   except Error as e:
      print(e)
   finally:
      if con:
         print('All tables in database up-to-date. Initialisation complete.')
         con.close()
      