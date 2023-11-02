import numpy as np
import pandas as pd
import sqlite3
import sys

__author__ = "Carolin Brune"

############################################################################
# variables
############################################################################

ALLOWED_DATABASE_LINKS = ['BiGG', 'MetaNetX', 'SEED', 'ChEBI', 'KEGG']
REQUIRED_SUBSTANCE_ATTRIBUTES = ['name', 'formula', 'flux', 'source']

############################################################################
# classes
############################################################################


class Medium:

    def __init__(self, name:str, substance_table:pd.DataFrame, description=None, doi=None):
        
        self.name = name
        self.description = description
        self.substance_table = substance_table
        self.doi = doi

    # -----------------------------------
    # TODO:
    # add functionalities from SPECIMEN ?
    # ------------------------------------

    def format_substance_table(self, format:str) -> pd.DataFrame:

        formatted_table = self.substance_table.copy()

        match format:
            # for enabling usage of this new medium class with the old growth functions
            case 'growth_old':
                formatted_table = formatted_table[['BiGG', 'flux']]
                formatted_table['mediumname'] = self.name

            # raise error for unknow input
            case _:
                raise ValueError(f'Unknown format type for table: {format}')
            
        return formatted_table

############################################################################
# functions for loading from DB
############################################################################

def load_substance_table_from_db(mediumname: str, database:str, type:str) -> pd.DataFrame:
    
    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    match type:
        # use this for debugging
        case 'testing':
            result = cursor.execute("""SELECT substance.id, substance.name, substance.formula, medium2substance.flux 
                                    FROM medium, medium2substance, substance
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id
                                    """, (mediumname,)) 
            substance_table = result.fetchall()
        # create table for documentation
        case 'documentation':
            result = cursor.execute("""SELECT substance.name, substance.formula, medium2substance.flux , medium2substance.source, substance2db.db_id, substance2db.db_type
                                    FROM medium, medium2substance, substance, substance2db
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id AND substance2db.substance_id = substance.id
                                    """, (mediumname,)) 
            substance_table = result.fetchall()
        # create table with all information, standard for generating the Medium table
        case 'standard':
            result = cursor.execute("""SELECT substance.name, substance.formula, medium2substance.flux , medium2substance.source, substance2db.db_id, substance2db.db_type
                                    FROM medium, medium2substance, substance, substance2db
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id AND substance2db.substance_id = substance.id
                                    """, (mediumname,)) 
            substance_table = result.fetchall()
        # default: throw error
        case _:
            connection.close()
            sys.exit("Unknown type for loading the substance table.")

    # close connection
    connection.close()

    return substance_table


def load_medium_from_db(name:str, database:str, type='standard') -> Medium:
    
    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # get description
    result = cursor.execute("SELECT medium.description FROM medium WHERE medium.name = ?",(name,))
    if result:
        description = result.fetchone()[0] # each name should be unique
    else:
        print('ERROR: Unknown medium name', file=sys.stderr)
        sys.exit(1)

    # get DOI(s)
    result = cursor.execute("SELECT medium.reference FROM medium WHERE medium.name = ?",(name,))
    doi = result.fetchone()[0]

    # close connection to DB
    connection.close()
    
    # get substance table
    substance = load_substance_table_from_db(name, database, type)

    return Medium(name=name, substance_table=substance, description=description, doi=doi)


############################################################################
# functions for reading media from extern
############################################################################

def read_substances_from_file(path: str) -> pd.DataFrame: 
    """...

    Format of the TSV:
    name | formula | flux | source | X | X | ...

    X: placeholder for database names (columns filled with corresponding IDs of the substances)
    X = see ALLOWED_DATABASE_LINKS

    Returns:
        _type_: _description_
    """

    # read in the table
    substance_table = pd.read_csv(path, sep='\t')

    # validate the table
    substance_cols = substance_table.columns
    for req in REQUIRED_SUBSTANCE_ATTRIBUTES:
        if req in substance_cols:
            continue
        else:
            substance_table[req] = None

    # remove unwanted database links
    all_allowed = REQUIRED_SUBSTANCE_ATTRIBUTES + ALLOWED_DATABASE_LINKS
    for_removal = [_ for _ in substance_cols if _ not in all_allowed] 
    substance_table.drop(for_removal, axis=1, inplace=True)

    # convert all NaNs into None
    substance_table = substance_table.replace({np.nan: None})

    return substance_table


def read_external_medium(how:str) -> Medium:

    match how:
        # interactive mode using the console
        case 'console':

            # prompt for name
            name = input('Enter a medium name:\n')

            # prompt for description
            description = input('Enter a short description for the medium:\n')

            # promt for description 
            reference = input('Enter the reference linke (DOI) for the medium:\n')

            # prompt for path to substance table
            substance_path = input('Enter a path to a table of substances (TSV):\n')
            # tsv format
            substances = read_substances_from_file(substance_path)

            # construct medium
            new_medium = Medium(name=name, description=description, substance_table=substances, doi=reference)

            return new_medium
        
        # .......................................................
        # TODO
        #    construct a case to read from a file or a list or so
        #    to read in multiple media from a list
        # .......................................................

        # unknown case, raise error
        case _:
            raise ValueError(f'Unknown input for parameter how: {how}')


############################################################################
# functions for adding entries to DB
############################################################################

def get_last_idx_table(tablename, connection, cursor):

    match tablename:
        case 'medium':
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM medium")
        case 'medium2substance':
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM medium2substance")
        case 'substance':
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM substance")
        case 'substance2db':
            last_rowid_res = cursor.execute("SELECT MAX(rowid) FROM substance2db")
        case _:
            sys.exit(f'Unknown table name for database: {tablename}')

    last_rowid_fetched = last_rowid_res.fetchone()
    last_rowid = last_rowid_fetched[0]

    return last_rowid


def enter_substance_row(row: pd.Series, connection, cursor):

    # check if a perfect match exists in database
    exact_match_res = cursor.execute("SELECT * FROM substance WHERE substance.formula = ? AND substance.name = ? ",(row['formula'],row['name']))
    exact_match = exact_match_res.fetchone()

    if exact_match:

        # get substance id
        substance_id = exact_match[0]

    # if not promt for similar or create new substance
    else:    

        candidates = cursor.execute("SELECT * FROM substance WHERE substance.formula = ? OR substance.name LIKE ? ",(row['formula'],row['name']))
        candidate_list = candidates.fetchall()

        # if candidates were found, let user chose if one of those should be used
        
        substance_id = None
        if len(candidate_list) > 0:

            print(f'Following similar entries have been found for the substance {row["name"]}.')
            for i in candidate_list:
                print(i)
            res = input('If one matches your entry, please enter the ID or enter skip/s: \n')
            while True:
                # skip
                if res in ['skip','s']:
                    break
                # match found
                elif res.isnumeric() and int(res) in [_[0] for _ in candidate_list]:
                    substance_id = int(res)
                    break
                # unknown input, try again
                else:
                    res = input('Unknow input. Please try again:\n')

        # if no matching entry has been found, insert new entry into DB
        if not substance_id:
            cursor.execute("INSERT INTO substance VALUES(?,?,?)",(None,row['name'],row['formula'],))
            connection.commit()
            # retrieve the newly added ID
            res = cursor.execute("SELECT last_insert_rowid() FROM substance")
            substance_id = res.fetchone()[0]
        
    return substance_id


def enter_m2s_row(row: pd.Series, medium_id, connection, cursor):

    cursor.execute('INSERT INTO medium2substance VALUES(?,?,?,?)',(medium_id,row['substance_id'],row['flux'],row['source'],))
    connection.commit()
    

def enter_s2db_row(row: pd.Series, db_type: str, connection, cursor):

    # only need to enter if ID exists
    if row[db_type]:
        # check if mapping already exists
        entry_check_res = cursor.execute("SELECT 1 FROM substance2db WHERE substance2db.substance_id = ? AND db_id = ?", (row['substance_id'],row[db_type],))
        entry_check = entry_check_res.fetchone()
        if not entry_check:
            cursor.execute('INSERT INTO substance2db VALUES(?,?,?)',(row['substance_id'],row[db_type],db_type,))
            connection.commit()


def enter_medium_into_db(database: str, medium: Medium):

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # collect starting IDs
    sid_medium = get_last_idx_table('medium', connection, cursor)
    sid_medium2substance = get_last_idx_table('medium2substance', connection, cursor)
    sid_substance = get_last_idx_table('substance', connection, cursor)
    sid_substance2db = get_last_idx_table('substance2db', connection, cursor)

    # try to add a new medium to the database
    try:

        # name
        #   check, if name already in DB and prompt, if that problem occurs
        #   enter into database, if no problem remains 
        name_check_res = cursor.execute("SELECT 1 FROM medium WHERE medium.name = ?", (medium.name,))
        name_check = name_check_res.fetchone()
        if not name_check:
            # medium name is unique 
            # add medium to DB
            cursor.execute("INSERT INTO medium VALUES(?,?,?,?)",(None,medium.name,medium.description,medium.doi,))
        else:
            # medium name already in DB 
            check = input('The name {medium.name} is already in the database.\nDo you want to set a new name? (yes/no): ')
            if check in ['y','yes']:
                # set a new name
                while True:
                    new_name = input('Please enter a medium name:\n')
                    name_check_2 = cursor.execute("SELECT 1 FROM medium WHERE medium.name = ?", (new_name,))
                    if name_check_2.fetchone():
                        print('This name already exists in the database.')
                    else:
                        medium.name = new_name
                        break
            elif check in ['n','no']:
                # end program when no new name is set
                print('No new name chosen. Ending the program.')
                sys.exit()
            else:
                # Abort program at unknown input
                sys.exit('Unknown input. Aborting.')
        
        connection.commit()
        res = cursor.execute("SELECT last_insert_rowid() FROM medium")
        medium_id = res.fetchone()[0]

        # substances 
        medium.substance_table['substance_id'] = medium.substance_table.apply(enter_substance_row, axis=1, args=(connection, cursor))

        # medium2substances
        medium.substance_table.apply(enter_m2s_row, axis=1, args=(medium_id, connection, cursor))

        # substance2db
        avail_data = [_ for _ in medium.substance_table.columns if _ in ALLOWED_DATABASE_LINKS]
        for current_db in avail_data:
            medium.substance_table.apply(enter_s2db_row, axis=1, args=(current_db,connection,cursor))

        print(f'Medium {medium.name} with ID {medium_id} has been added to the database.')

    # in case of problems, interupt and revert changes
    except:

        print(f'During execution, the following error occured: \n{sys.exc_info()}')
        print('Reverting changes to database...')

        cursor.execute("DELETE FROM medium WHERE rowid > ?",(sid_medium,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?",(sid_medium2substance,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?",(sid_substance,))
        cursor.execute("DELETE FROM medium WHERE rowid > ?",(sid_substance2db,))

        print('Done')

    # in any case close connection to database
    finally:
        connection.close()


############################################################################
# possible entry points
############################################################################

# entry point for entering a medium using the command line
# TODO since database is part of package, direct accessing possible
def add_medium(database:str):
    
    # get external medium
    medium = read_external_medium('console')

    # add to database
    enter_medium_into_db(database, medium)