import io
import cobra
import copy
import numpy as np
import pandas as pd
from pathlib import Path
import sqlite3
import sys
import warnings
from typing import Literal, Union
import random
import string
from sqlite_dump import iterdump

__author__ = "Carolin Brune"

############################################################################
# variables
############################################################################

PATH_TO_DB = Path(Path(__file__).parent.resolve(), 'data.db')

ALLOWED_DATABASE_LINKS = ['BiGG', 'MetaNetX', 'SEED', 'VMH', 'ChEBI', 'KEGG']
REQUIRED_SUBSTANCE_ATTRIBUTES = ['name', 'formula', 'flux', 'source']

# subset addable to media
# -----------------------

CAS_DB_NAMES = ['L-Alanine', 'L-Arginine', 'L-Aspartate', 'L-Glutamate', 'Glycine', 'L-Histidine', 'L-Isoleucine', 'L-Leucine', 
                'L-Lysine', 'L-Methionine', 
                'L-Phenylalanine', 'L-Proline', 'L-Serine','L-Threonine','L-Tryptophan', 'L-Tyrosine', 'L-Valine', 'L-Cystine']
PROTEINOGENIC_AA_DB_NAMES = ['L-Alanine', 'L-Arginine', 'L-Asparagine', 'L-Aspartate', 'L-Cysteine', 
                             'L-Glutamine', 'L-Glutamate', 'Glycine', 'L-Histidine', 'L-Isoleucine', 
                             'L-Leucine', 'L-Lysine', 'L-Methionine', 'L-Phenylalanine', 'L-Proline', 
                             'L-Serine', 'L-Threonine', 'L-Tryptophan', 'L-Tyrosine', 'L-Valine']
SUBSET_MEDIA_MAPPING = {'casamino':CAS_DB_NAMES, 
                        'aa':PROTEINOGENIC_AA_DB_NAMES}


############################################################################
# classes
############################################################################


class Medium:
    """Class describing a medium.

    Attributes:
        name (str): The name or abbreviation of the medium.
        substance_table (pd.DataFrame): A table containing information about the medium in silico compounds. Long format.
        description (str, optional): Short description of the medium. Defaults to None.
        doi (str): Reference(s) to the original publication of the medium. Defaults to None.
    """

    def __init__(self, name:str, substance_table:pd.DataFrame=pd.DataFrame(columns=['name','formula','flux','source','db_id','db_type']), description:str=None, doi:str=None):
        """Initialise a Medium object.

        Args:
            name (str): The name or abbreviation of the medium.
            substance_table (pd.DataFrame, optional):  A table containing information about the medium in silico compounds. Long format.
                Defaults to an empty table with the columns ['name','formula','flux','source','db_id','db_type'].
            description (str, optional): Short description of the medium.. Defaults to None.
            doi (str, optional): Reference(s) to the original publication of the medium.. Defaults to None.
        """
        
        self.name = name
        self.description = description
        self.substance_table = substance_table
        self.doi = doi

    # -----------------------------------
    # TODO:
    # add functionalities from SPECIMEN ?

    # possible @TODO
    # ..............
    # has compound
    # remove compound
    # set source of
    # ------------------------------------

    def add_substance_from_db(self, name:str, flux:float=10.0):
        """Add a substance from the database to the medium.

        Args:
            name (str): Name of the substance. Should be part of the database substance.name column.
            flux (float, optional): Sets the flux value of the new substance. Defaults to 10.0.
        """
        # build connection to DB
        connection = sqlite3.connect(PATH_TO_DB)
        cursor = connection.cursor()
        # fetch substance
        result = cursor.execute(""" SELECT substance.name, substance.formula, substance2db.db_id, substance2db.db_type 
                                    FROM substance, substance2db 
                                    WHERE substance.name = ? AND substance.id = substance2db.substance_id""", (name,))
        substance = result.fetchall()

        # check if fetch was successful
        if len(substance) == 0:
            warnings.warn(f'Could not fetch substance {name} from DB.')
            return
        
        # add substance to table
        substance_table = pd.DataFrame(substance, columns=['name','formula','db_id','db_type'])
        substance_table.insert(2,'flux',flux)
        substance_table.insert(3,'source',None)
        self.substance_table = pd.concat([self.substance_table, substance_table], ignore_index=True)


    def get_source(self, element:str) -> list[str]:
        """Get the source of a given element for the medium.

        Search for the given element (elemental symbol e.g. O), excluding pattern matches that are 
        followed by other lower-case letters and returm them as a list of sources for the given element. 

        Args:
            element (str): The symbol of the element to search the sources of

        Returns:
            list[str]: The list of the names of the sources (no duplicates).
        """

        return list(set(self.substance_table[self.substance_table['formula'].str.contains(element + '(?![a-z])', case=True, regex=True)]['name']))


    def is_aerobic(self) -> bool:
        """Check if the medium contains O2 / dioxygen.

        Returns:
            bool: Results of the test, True if pure oxygen is in the medium.
        """

        test = self.substance_table[['name','formula']].loc[(self.substance_table['name'] == 'Dioxygen [O2]') & (self.substance_table['formula'] == 'O2')].any().all()
        return test
    
    
    def make_aerobic(self, flux:float=None):
        """If the medium is curretly anaerobic, add oxygen to the medium to make it aerobic.
        
        Args:
            flux(float,optional): The flux value for the oxygen to be added. Defaults to None.
        """
        if not self.is_aerobic():
            # build connection to DB
            connection = sqlite3.connect(PATH_TO_DB)
            cursor = connection.cursor()
            # fetch oxygen
            result = cursor.execute(""" SELECT substance.name, substance.formula, substance2db.db_id, substance2db.db_type 
                                    FROM substance, substance2db 
                                    WHERE substance.formula = 'O2' AND substance.name = 'Dioxygen [O2]' AND substance.id = substance2db.substance_id""")
            oxygen = result.fetchall()
            # add oxygen to substance table
            oxygen_table = pd.DataFrame(oxygen, columns=['name','formula','db_id','db_type'])
            oxygen_table.insert(2,'flux',flux)
            oxygen_table.insert(3,'source',None)
            self.substance_table = pd.concat([self.substance_table, oxygen_table], ignore_index=True)


    def make_anaerobic(self):
        """If the medium is currently aerobic, deletes the oxygen from it to make it anaerobic.
        """
        if self.is_aerobic():
            # remove dioxygen // O2 from the substance table
            self.substance_table.drop(self.substance_table[(self.substance_table['name']=='Dioxygen [O2]') & (self.substance_table['formula']=='O2')].index, inplace=True)

    # @IDEA: provide dict with names and fluxes to set some of them separatly if user wishes to do so - or separate function?
    def set_default_flux(self,flux:float =10.0, replace:bool =False, double_o2:bool =True):
        """Set a default flux for the model.

        Args:
            flux (float, optional): Default flux for the medium. Defaults to 10.0.
            replace (bool, optional): 
            double_o2 (bool, optional): Tag to double the flux for oxygen only. Works only with replace=True. Defaults to True.
        """
        # replace fluxes
        if replace:
            self.substance_table['flux'] = flux
            if self.is_aerobic() and double_o2:
                self.substance_table.loc[self.substance_table['name']=='Dioxygen [O2]','flux'] = 2*flux
        # keep already set fluxes
        else:
            self.substance_table['flux'] = self.substance_table['flux'].fillna(flux)
            if self.is_aerobic() and double_o2 and self.substance_table.loc[self.substance_table['name']=='Dioxygen [O2]','flux'].any():
                self.substance_table.loc[self.substance_table['name']=='Dioxygen [O2]','flux'] = 2*flux
        
    def set_oxygen_percentage(self, perc:float=1.0):
        """Set oxygen percentage of the medium.

        Args:
            perc (float, optional): Percentage of oxygen. Defaults to 1.0 (= 100%)
        """

        if self.is_aerobic():
            self.substance_table.loc[self.substance_table['name']=='Dioxygen [O2]','flux'] = self.substance_table.loc[self.substance_table['name']=='Dioxygen [O2]','flux'] * perc
        else:
            warnings.warn(f'WARNING: no oxygen detected in medium {self.name}, cannot set oxygen percentage.')


    def combine(self,other:'Medium') -> 'Medium':
        """Combine two media into a new one.

        Args:
            other (Medium): A second medium the first one is to be combined with.

        Returns:
            Medium: The combined medium, without fluxes and sources.
        """

        combined = copy.deepcopy(self)
        combined.name = self.name + '+' + other.name
        combined.description = f'Combined medium contructed from {self.name} and {other.name}.'
        combined.substance_table = pd.concat([combined.substance_table, other.substance_table], ignore_index=True)
        combined.doi = str(self.doi) + ', ' + str(other.doi)

        # remove fluxes and sources, as they are no longer a fit
        combined.substance_table['flux'] = None
        combined.substance_table['source'] = None

        # remove duplicate rows
        combined.substance_table.drop_duplicates(inplace=True,ignore_index=True)

        return combined
    
    
    def __add__(self,other:'Medium') -> 'Medium':
        return self.combine(other)

    def add_subset(self, type: Literal['aa', 'casamino']) -> 'Medium':
        """Add a subset of substances to the medium, returning a newly generated one.
        Available subset are the following:
        - aa: the 20n proteinogenic amino acids
        - casamino: casamino acid, based on USBiological Life Sciences

        Args:
            type(str): The type of subset to be added. Choices are 'aa' and 'casamino'.

        Returns:
            Medium: A new medium that is the combination of the set subset and the old one.
        """

        # get substance names based on type
        subset_names = SUBSET_MEDIA_MAPPING[type]

        # open connection to database
        connection = sqlite3.connect(PATH_TO_DB)
        cursor = connection.cursor()

        # extract amino acids from DB
        substances = []
        for name in subset_names:
            result = cursor.execute("""SELECT 1 FROM substance WHERE substance.name = ?""",(name,))
            check = result.fetchone()
            if check:
                result = cursor.execute("""SELECT substance.name, substance.formula, substance2db.db_id, substance2db.db_type
                                        FROM substance, substance2db
                                        WHERE substance.name = ? AND substance.id = substance2db.substance_id""",(name,))
                subs = result.fetchone()
                substances.append(subs)
            else:
                warnings.warn(f'Could not find substance in DB: {name}')

        # reformat into Medium object
        substances = pd.DataFrame(substances, columns=['name','formula','db_id','db_type'])
        substances.insert(2,'flux',None)
        substances.insert(3,'source',None)

        aa_medium = Medium(name=type, substance_table=substances, description=type)

        # combine with current one 
        return self + aa_medium 


    # functions for export table
    # --------------------------

    # @TODO
    def format_substance_table(self, format:str) -> pd.DataFrame:
        """Produce a reformatted version of the substance table for different purposes.

        Possible formats
        - 'growth_old': for working with the old version of the growth module.

        Args:
            format (str): Specifies the output format type.

        Raises:
            ValueError: Unknown format type for table.

        Returns:
            pd.DataFrame: The reformatted substance table (copy).
        """

        match format:
            # for enabling usage of this new medium class with the old growth functions
            case 'growth_old':

                formatted_table = self.substance_table.loc[self.substance_table['db_type']=='BiGG'][['name','flux','db_id']]
                formatted_table.rename(columns={'db_id':'BiGG'})

                # ..............................................................
                # TODO / WARNING
                # currently only substances WITH a BiGG ID are kept in this step
                # ..............................................................
                
                formatted_table['mediumname'] = self.name

            # raise error for unknow input
            case _:
                raise ValueError(f'Unknown format type for table: {format}')
            
        return formatted_table
    

    def produce_medium_docs_table(self, folder: str = './', max_width: int = 80) -> str:
        """Produces a rst-file containing reStructuredText for the substance table for documentation.

        Args:
            folder (str, optional): Path to folder/directory to save the rst-file to. Defaults to './'.
            max_width (int, optional): Maximal table width of the rst-table. Defaults to 80.
        """

        def calculate_column_widths_docs(header: list, max_width: int, flux_width = 8) -> str:
            """Helper function for :py:func:`produce_medium_docs_table`. 
            Calculates the columns widths based on the header lengths and maximal widths.

            Args:
                header (list): List of column names of the table.
                max_width (int): Maximal widths of the output table.
                flux_width (int, optional): Widths of the flux column, if 'flux' is part of 'header'. Defaults to 8.
 
            Returns:
                str: The columns widths as a string, separated by a whitespace.
            """

            if len(header) == 2:
                partition = max_width // 3
                return f"{str(max_width-partition)} {partition}"
            else:
                partition = (max_width - flux_width) // 2
                return f"{str(max_width-flux_width-partition)} {flux_width} {partition}"

        def produce_medium_docs_table_row(row: pd.Series, file: io.TextIOWrapper):
            """Helper function for :py:func:`produce_medium_docs_table`.
            Tranforms each row of the substance table into a row of the rst-file.

            Args:
                row (pd.Series): The row of the Medium.substance_table.
                file (io.TextIOWrapper): The connection to the file to write the rows into.
            """

            list = row.to_list()
            file.write(f"  * - {list[0]}\n")
            for l in list[1:]:
                file.write(f"    - {l}\n")

        with open(folder + f'{self.name}.rst', 'w') as f:

            print(type(f))

            # slim table to columns of interest for documentation
            m_subs = self.substance_table[['name','flux','source']]

            if all(m_subs['flux'].values == None):
                m_subs = m_subs.drop('flux', axis=1)

            header = m_subs.columns

            widths = calculate_column_widths_docs(header, max_width)
            
            title = f'{self.description} ({self.name})'
            
            # Produce header/title of HTML page
            f.write(
                f"{title}\n"
                f"{'^' * len(title)}\n\n"
            )

            # produce descriptor
            f.write(f".. list-table:: {self.name} composition\n"
                    f"  :name: {self.name.lower()}_comp\n"
                    "  :align: center\n"
                    f"  :widths: {widths}\n"
                    "  :header-rows: 1\n"
                    "  :class: no-scrollbar-table\n\n")
            
            # produce header
            f.write(f"  * - {header[0]}\n")
            for l in header[1:]:
                f.write(f"    - {l}\n")

            # produce table body
            m_subs.apply(produce_medium_docs_table_row, file=f, axis=1)

            f.close()


    def export_to_file(self, type:str='tsv',dir:str='./', max_widths: int = 80):
        """Export medium, especially substance table.

        Args:
            type (str, optional): Type of file to export to. Defaults to 'tsv'. Further choices are 'csv', 'docs', 'rst'.
            dir (str, optional): Path to the directory to write the file to. Defaults to './'.
            max_widths (int, optional): Maximal table width for the documentation table (). Only viable for 'rst' and 'docs'.
                Defaults to 80.

        Raises:
            ValueError: Unknown export type if type not in ['tsv','csv','docs','rst']
        """

        path_no_ext = dir + self.name + '_substances'

        match type:
            case 'tsv':
                self.substance_table.to_csv(path_no_ext + '.tsv', sep='\t', index=False)
            case 'csv':
                self.substance_table.to_csv(path_no_ext + '.csv', sep=';', index=False)
            case 'docs' | 'rst':
                self.produce_medium_docs_table(folder = dir, max_width = max_widths)
            case _:
                raise ValueError('Unknown export type: {type}')
            

    # functions for conversion
    # ------------------------
    # @TODO
    #    extend namespace options
    #        maybe a warning or something similar, if compounds do not have an ID corresponding to the namespace
    # @TEST
    #    name option for namespace
    # @DEV : extent literals in all related functions after extension of namespace options
    def export_to_cobra(self, namespace:Literal['Name', 'BiGG']='BiGG', default_flux:float=10.0, replace:bool=False, double_o2:bool=True) -> dict[str,float]:
        
        self.set_default_flux(default_flux, replace=replace, double_o2=double_o2)

        match namespace:
            case 'Name':

                names = self.substance_table[['name','flux']]
                names['EX'] = biggs['db_id'] + ' exchange'
                cobra_medium = pd.Series(names.flux.values,index=names.EX).to_dict()       

            case 'BiGG':
                
                biggs = self.substance_table[self.substance_table['db_type']=='BiGG'][['name','db_id','flux']]
                biggs['db_id_EX'] = 'EX_' + biggs['db_id'] + '_e'
                cobra_medium = pd.Series(biggs.flux.values,index=biggs.db_id_EX).to_dict()       

                # .....................................
                # TODO
                #    what about those without BiGG IDs?
                # .....................................       

            case _:
                raise ValueError(f'Unknown namespace: {namespace}')
        
        return cobra_medium  

############################################################################
# functions for loading from DB
############################################################################

def load_substance_table_from_db(mediumname: str, database:str, type='standard') -> pd.DataFrame:
    """Load a substance table from a database.

    Currently available types:
    - 'testing': for debugging
    - 'documentation': for downloading a table for the docs
    - 'standard': The standard format containing all information in long format.

    Note: 'documentation' currently object to change

    Args:
        name (str): The name (or identifier) of the medium.
        database (str): Path to the database.
        type (str, optional): How to load the table. Defaults to 'standard'.

    Raises:
        ValueError: Unknown type for loading the substance table.

    Returns:
        pd.DataFrame: The substance table in the specified type retrieved from the database.
    """
    
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
            substance_table = pd.DataFrame(substance_table, columns=['id','name','formula','flux'])

        # create table for documentation
        case 'documentation':
            result = cursor.execute("""SELECT substance.name, substance.formula, medium2substance.flux , medium2substance.source, substance2db.db_id, substance2db.db_type
                                    FROM medium, medium2substance, substance, substance2db
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id AND substance2db.substance_id = substance.id
                                    """, (mediumname,)) 
            substance_table = result.fetchall()
            substance_table = pd.DataFrame(substance_table, columns=['name','formula','flux','source','db_id','db_type'])

        # create table with all information, standard for generating the Medium table
        case 'standard':
            result = cursor.execute("""SELECT substance.name, substance.formula, medium2substance.flux , medium2substance.source, substance2db.db_id, substance2db.db_type
                                    FROM medium, medium2substance, substance, substance2db
                                    WHERE medium.name = ? AND medium.id = medium2substance.medium_id AND medium2substance.substance_id = substance.id AND substance2db.substance_id = substance.id
                                    """, (mediumname,)) 
            substance_table = result.fetchall()
            substance_table = pd.DataFrame(substance_table, columns=['name','formula','flux','source','db_id','db_type'])

        # default: throw error
        case _:
            connection.close()
            raise ValueError(f"Unknown type for loading the substance table: {type}")

    # close connection
    connection.close()

    return substance_table


def load_medium_from_db(name:str, database:str=PATH_TO_DB, type:str='standard') -> Medium:
    """Load a medium from a database.

    Args:
        name (str): The name (or identifier) of the medium.
        database (str, optional): Path to the database. Defaults to the in-built database.
        type (str, optional): How to load the medium. Defaults to 'standard'.

    Raises:
        ValueError: Unknown medium name.

    Returns:
        Medium: The medium retrieved from the database.
    """
    
    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # get description
    result = cursor.execute("SELECT medium.description FROM medium WHERE medium.name = ?",(name,))
    if result:
        description = result.fetchone()[0] # each name should be unique
    else:
        raise ValueError(f'Unknown medium name: {name}')

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
    """Read in a TSV with substance information into a table.

    Format of the TSV:
    name | formula | flux | source | X | X | ...

    X: placeholder for database names (columns filled with corresponding IDs of the substances)
    X = see ALLOWED_DATABASE_LINKS

    Returns:
        pd.DataFrame: The table of substance information read from the file
    """

    # read in the table
    substance_table = pd.read_csv(path, sep='\t', comment='#')

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

# @TEST file flag for how
def read_external_medium(how:str, **kwargs) -> Medium:
    """Read in an external medium. 

    Currently available options for how:
    - 'console': read in medium via the console
    - 'file': read in medium from a file, requires a 'path=str' argument to be passed.

    About the format (console, file):
    The substances have to be saved in a TSV file table (format see :py:func:`read_substances_from_file`).
    Further information for the 'file' option have to added as comments in the format: `# info_name: info`.
    Information should but do not need to contain name, description and reference.

    Args:
        how (str): How (or from where) the medium should be read in.
            Available options are given above.

    Raises:
        ValueError: Unknown description for how.

    Returns:
        Medium: The read-in medium.
    """

    match how:
        # interactive mode using the console
        case 'console':

            # prompt for name
            name = input('Enter a medium name:\n')

            # prompt for description
            description = input('Enter a short description for the medium:\n')

            # promt for description 
            reference = input('Enter the reference link (DOI) for the medium:\n')

            # prompt for path to substance table
            substance_path = input('Enter a path to a table of substances (TSV):\n')
            # tsv format
            substances = read_substances_from_file(substance_path)

            # construct medium
            new_medium = Medium(name=name, description=description, substance_table=substances, doi=reference)

            return new_medium
        
        # read from a file 
        case 'file':

            # get info from comments
            infos = {}
            with open(kwargs['path'], 'r') as f:
                for line in f:
                    # parse the comment lines from the file
                    if line.startswith('#'):
                        l = line.split(':')
                        k = l[0]
                        k = k[1:].strip()
                        v = "".join(l[1:])
                        v = v.strip()
                        infos[k] = v
                    else:
                        break
            # retrieve important infos and supplement missing
            if 'name' in infos.keys():
                name = infos['name']
            else:
                warnings.warn('No name found for externally loaded medium. Setting random one.')
                name = 'noname_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
            
            if 'description' in infos.keys():
                description = infos['description']
            else:
                description = None
            
            if 'reference' in infos.keys():
                reference = infos['reference']
            else:
                reference = None

            # get substances
            substances = read_substances_from_file(kwargs['path'])

            # construct medium
            new_medium = Medium(name=name, description=description, substance_table=substances, doi=reference)

            return new_medium

        # unknown case, raise error
        case _:
            raise ValueError(f'Unknown input for parameter how: {how}')


def extract_medium_info_from_model_bigg(row, model:cobra.Model) -> pd.Series:
    """Helper function for :py:func:`read_from_cobra_model`. 
    Extracts more information about the medium.

    Args:
        row (pd.Series): A row of the datatable of :py:func:`read_from_cobra_model`.
        model (cobra.Model): The cobra Model

    Returns:
        pd.Series: _description_
    """

    db_id = row['db_id_EX'].replace('EX_','').replace('_e','')
    meta = model.metabolites.get_by_id(db_id+'_e')
    name = meta.name
    formula = meta.formula

    return pd.Series([name,formula,db_id])


def read_from_cobra_model(model: cobra.Model) -> Medium:
    """Read and import a medium from a cobra model into a Medium object.

    Args:
        model (cobra.Model): An open cobra Model.

    Returns:
        Medium: The imported medium.
    """

    # retrieve the medium from the model
    cobra_medium = model.medium.copy()
    substances = pd.DataFrame(cobra_medium.items(),columns=['db_id_EX','flux'])
    substances['db_type'] = 'BiGG'

    # retrieve additional information
    substances['source'] = model.id
    substances['name'] = None
    substances.apply(extract_medium_info_from_model_bigg, model=model, axis=1)

    # reformat table 
    substances.drop(columns=['db_id_EX'],inplace=True)
    sorted_cols = ['name','formula','flux','source','db_id','db_type']
    substances = substances.reindex(columns=sorted_cols)

    # construct the medium 
    name = f'medium_{model.id}'
    description = f'Medium imported from model {model.id}'
    imported_medium = Medium(name, substances, description)

    return imported_medium


############################################################################
# functions for adding entries to DB
############################################################################

def get_last_idx_table(tablename: str, connection: sqlite3.Connection, cursor: sqlite3.Cursor) -> int:
    """Helper function for :py:func:`refinegems.database.enter_medium_into_db`.
    Retrieves the last row id of a specified table of the database.

    Args:
        tablename (str): The name of the table to retrieve the last row id from
        connection (sqlite3.Connection): Connection to the database.
        cursor (sqlite3.Cursor): Cursor for the database.

    Returns:
        int: The last row ID of the specified table.
    """

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


def enter_substance_row(row: pd.Series, connection: sqlite3.Connection, cursor:sqlite3.Cursor) -> int:
    """Helper function for :py:func:`refinegems.database.enter_medium_into_db`.
    Enters a new entry in the medium2substance table.

    Args:
        row (pd.Series): A row of the pd.DataFrame of the :py:func:`refinegems.database.enter_medium_into_db` function.
        connection (sqlite3.Connection): Connection to the database.
        cursor (sqlite3.Cursor): Cursor for the database.

    Returns:
        int: The substance ID in the database of the substance.
    """

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


def enter_m2s_row(row: pd.Series, medium_id: int, connection: sqlite3.Connection, cursor: sqlite3.Cursor):
    """Helper function for :py:func:`refinegems.database.enter_medium_into_db`.
    Enters a new entry in the medium2substance table.

    Args:
        row (pd.Series): A row of the pd.DataFrame of the :py:func:`refinegems.database.enter_medium_into_db` function.
        medium_id (int): The row id of the medium.
        connection (sqlite3.Connection): Connection to the database.
        cursor (sqlite3.Cursor): Cursor for the database.
    """

    # check if entry already exists in database
    exact_match_res = cursor.execute("SELECT 1 FROM medium2substance WHERE medium2substance.medium_id = ? AND medium2substance.substance_id = ? ",(medium_id,row['substance_id']))
    exact_match = exact_match_res.fetchone()

    # else
    if not exact_match:
        cursor.execute('INSERT INTO medium2substance VALUES(?,?,?,?)',(medium_id,row['substance_id'],row['flux'],row['source'],))
        connection.commit()
    else:
        warnings.warn(f'Medium2substance connection {medium_id} - {row["substance_id"]} already exists, skipped second assignment.')
    

def enter_s2db_row(row: pd.Series, db_type: str, connection: sqlite3.Connection, cursor: sqlite3.Cursor):
    """Helper function for :py:func:`refinegems.database.enter_medium_into_db`.
    Enters a new entry in the substance2db table after checking if it has yet to be added.

    Args:
        row (pd.Series): A row of the pd.DataFrame of the :py:func:`refinegems.database.enter_medium_into_db` function.
        db_type (str): Type of database identifier to be added.
        connection (sqlite3.Connection): Connection to the database.
        cursor (sqlite3.Cursor): Cursor for the database.
    """

    # only need to enter if ID exists
    if row[db_type]:
        # check if mapping already exists
        entry_check_res = cursor.execute("SELECT 1 FROM substance2db WHERE substance2db.substance_id = ? AND db_id = ?", (row['substance_id'],row[db_type],))
        entry_check = entry_check_res.fetchone()
        if not entry_check:
            cursor.execute('INSERT INTO substance2db VALUES(?,?,?)',(row['substance_id'],row[db_type],db_type,))
            connection.commit()


def enter_medium_into_db(database: str, medium: Medium, update:bool=False):
    """Enter a new medium to an already existing database.

    Args:
        database (str): Path to the database.
        medium (Medium): A medium object to be added to the database.
        update (bool): Specifies if an existing medium should be updated or a new medium is added 
    """

    # build connection to DB
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    # collect starting IDs
    sid_medium = get_last_idx_table('medium', connection, cursor)
    sid_medium2substance = get_last_idx_table('medium2substance', connection, cursor)
    sid_substance = get_last_idx_table('substance', connection, cursor)
    sid_substance2db = get_last_idx_table('substance2db', connection, cursor)

    # try to add a new/updated medium to the database
    try:

        # name
        #   check, if name already in DB and prompt, if that problem occurs
        #   enter into database, if no problem remains 
        name_check_res = cursor.execute("SELECT 1 FROM medium WHERE medium.name = ?", (medium.name,))
        name_check = name_check_res.fetchone()
        if not (name_check and update): # Medium name is unique & no update is requested
            # medium name is unique 
            # add medium to DB
            cursor.execute("INSERT INTO medium VALUES(?,?,?,?)",(None,medium.name,medium.description,medium.doi,))
        elif name_check and not update: # Medium name already exists & no update is requested
            # medium name already in DB 
            check_new = input('The name {medium.name} is already in the database.\nDo you want to set a new name? (yes/no): ')
            if check_new in ['y','yes']:
                # set a new name
                while True:
                    new_name = input('Please enter a medium name:\n')
                    name_check_2 = cursor.execute("SELECT 1 FROM medium WHERE medium.name = ?", (new_name,))
                    if name_check_2.fetchone():
                        print('This name already exists in the database.')
                    else:
                        medium.name = new_name
                        break
            elif check_new in ['n','no']:
                # end program when no new name is set
                print('No new name chosen. Ending the program.')
                sys.exit()
            else:
                # Abort program at unknown input
                sys.exit('Unknown input. Aborting.')        
        elif not name_check and update: # Medium name is unique, BUT update was requested
            # Medium name not in DB 
            check_update = input('The name {medium.name} is not in the database yet.\nDo you want to add a new medium instead? (yes/no): ')
            if check_update in ['y','yes']:
                enter_medium_into_db(database, medium)
            elif check_update in ['n','no']:
                # end program when user did not use the update flag incorrectly
                print('No medium can be updated with the provided medium. Ending the program.')
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
# working with models
############################################################################

# @Testing
def medium_to_model(model:cobra.Model, medium:Medium, namespace:str='BiGG', 
                    default_flux:float=10.0, replace:bool=False, double_o2:bool=True, 
                    add:bool=True) -> Union[None, dict[str, float]]:
    """Add a medium to a COBRApy model. 

    Args:
        model (cobra.Model): A model loaded with COBRApy.
        medium (Medium): A refinegems Medium object.
        namespace (str, optional): String to set the namespace to use for the model IDs. Defaults to 'BiGG'.
        default_flux (float, optional): Set a default flux for NaN values or all. Defaults to 10.0.
        replace (bool, optional): Option to replace existing flux values with the default if set to True. Defaults to False.
        double_o2 (bool, optional): Double the oxygen amount in the medium. Defaults to True.
        add (bool, optional): If True, adds the medium to the model, else the exported medium is returned. Defaults to True.

    Returns:
        Union[None, dict[str, float]]: Either none or the exported medium.
    """
    
    # export medium to cobra
    exported_medium = medium.export_to_cobra(namespace=namespace, default_flux=default_flux, replace=replace, double_o2=double_o2)

    # remove exchanges that do not exist in model
    model_exchanges = [_.id for _ in model.exchanges]
    exported_medium = {k:v for k,v in exported_medium.items() if k in model_exchanges}

    if add:
        # add to model
        model.medium = exported_medium
        return 
    else:
        return exported_medium


############################################################################
# possible entry points
############################################################################

# Function to extract SQL schema with updated SBO/media tables
def updated_db_to_schema():
    """Extracts the SQL schema from the database data.db & Transfers it into an SQL file
    """
    # Not needed to be included in Schema
    NOT_TO_SCHEMA = [
        'BEGIN TRANSACTION;', 'COMMIT;', 
        'bigg_to_sbo', 'ec_to_sbo',
        'bigg_metabolites', 'bigg_reactions', 
        'modelseed_compounds'
        ]
    counter = 0 # To count rows in newly generated file
    
    conn = sqlite3.connect(PATH_TO_DB)
    with open('./updated_sbo_media_db.sql', 'w') as file:
        for line in iterdump(conn):
            if not (any(map(lambda x: x in line, NOT_TO_SCHEMA))):
                if 'CREATE TABLE' in line and counter != 0:
                    file.write(f'\n\n{line}\n')
                else: file.write(f'{line}\n')
                counter += 1
    conn.close()


# entry point for entering a medium using the command line
# TODO since database is part of package, direct accessing possible
def add_medium(database:str):
    
    # get external medium
    medium = read_external_medium('console')

    # add to database
    enter_medium_into_db(database, medium)
    
    # Generate updated SQl schema
    updated_db_to_schema()