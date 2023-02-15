#!/usr/bin/env python
import io
import sqlite3
import requests
import pandas as pd
from sqlite3 import Error
from typing import Literal
from os import path

__author__ = 'Gwendolyn O. Gusak'


PATH_TO_DB_DATA = path.join(path.abspath(path.dirname(path.dirname(__file__))), 'data/databases')
PATH_TO_DB = path.join(PATH_TO_DB_DATA, 'data.db')
VERSION_FILE = path.join(PATH_TO_DB_DATA, 'current_bigg_db_version.txt')
VERSION_URL = 'http://bigg.ucsd.edu/api/v2/database_version'


def is_valid_database(db_cursor: sqlite3.Cursor) -> Literal['BIGG', 'SBO', 'COMPLETE', 'EMPTY']:
   """
   Verifies if database has:
      - 2 tables with names 'bigg_metabolites' & 'bigg_reactions'
      - 2 tables with names 'bigg_to_sbo' & 'ec_to_sbo'
   
   Args:
      db_cursor: sqlite3.Cursor object of a database
 
   Returns:
      Literal:
         'BIGG': Only BiGG tables are in data.db
         'SBO': Only SBO tables are in data.db
         'COMPLETE': All tables are in data.db
         'EMPTY': data.db is either empty or incorrect
   """
   # Fetches the table names as string tuples from the connected database
   db_cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
   tables = [string[0] for string in db_cursor.fetchall()]
   
   bigg_tables_contained = 'bigg_metabolites' in tables and 'bigg_reactions' in tables
   sbo_tables_contained = 'bigg_to_sbo' in tables and 'ec_to_sbo' in tables
   
   if bigg_tables_contained and len(tables) == 2: return 'BIGG'
   elif sbo_tables_contained and len(tables) == 2: return 'SBO'
   elif bigg_tables_contained and sbo_tables_contained and len(tables) == 4: return 'COMPLETE'
   else: return 'EMPTY'


def create_SBO_database(db_cursor: sqlite3.Cursor):
   """Creates the SBO annotation database with 2 tables ('bigg_to_sbo' & 'ec_to_sbo')
      from file './data/databases/sbo_database.sql'
   """
   with open(path.join(PATH_TO_DB_DATA, 'sbo_database.sql')) as schema:
      db_cursor.executescript(schema.read())


def update_bigg_db(latest_version: str, db_connection: sqlite3.Connection):
   """
   """
   print('Adding BiGG tables...')
   
   # Store currently used version
   with open(VERSION_FILE, 'w') as file:
      file.write(latest_version)
   
   # Remove BiGG tables if exist
   db_cursor = db_connection.cursor()
   db_cursor.execute('DROP TABLE IF EXISTS bigg_metabolites')
   db_cursor.execute('DROP TABLE IF EXISTS bigg_reactions')
   
   # Create BiGG metabolites table
   bigg_models_metabs_url = 'http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt'
   bigg_models_metabs = requests.get(bigg_models_metabs_url).text
   bigg_models_metabs_df = pd.read_csv(io.StringIO(bigg_models_metabs), sep='\t')
   bigg_models_metabs_df.to_sql('bigg_metabolites', db_connection, index=False)

   # Create BiGG reactions table
   bigg_models_reacs_url = 'http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
   bigg_models_reacs = requests.get(bigg_models_reacs_url).text
   bigg_models_reacs_df = pd.read_csv(io.StringIO(bigg_models_reacs), sep='\t')
   bigg_models_reacs_df.to_sql('bigg_reactions', db_connection, index=False)
   

def get_latest_bigg_databases(db_connection: sqlite3.Connection, is_missing: bool=True):
   """
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
    
         
def initialise_database():
   
   # Initialise empty connection
   con = None

   print('Initialising database...')
   
   # Try to open connection & get cursor
   try:
      con = sqlite3.connect(PATH_TO_DB)
      cursor = con.cursor()
      
      print('Verifying database...')
      validity_code = is_valid_database(cursor)
      
      if validity_code == 'BIGG':
         print('Only BiGG tables contained in database. Adding SBO tables...')
         create_SBO_database(cursor)
      
      elif validity_code == 'SBO':
         print('Only SBO tables contained in database. Adding BiGG tables...')
         get_latest_bigg_databases(con)
         
      elif validity_code == 'EMPTY':
         print('Incorrect or empty database. Initialise database with required tables...')
         print('Adding SBO tables...')
         create_SBO_database(cursor)
         get_latest_bigg_databases(con)
         
      elif validity_code == 'COMPLETE':
         print('Verifying if BiGG tables are up-to-date...')
         get_latest_bigg_databases(con, False)
      
   except Error as e:
      print(e)
   finally:
      if con:
         print('All tables in database up-to-date. Initialisation complete.')
         con.close()
      