#!/usr/bin/env python
""" Provides functions to load and write models, media definitions and the manual annotation table

Depending on the application the model needs to be loaded with cobra (memote) or with libSBML (activation of groups). 
The media definitions are denoted in a csv within the data folder of this repository, thus the functions will only work if the user clones the repository. 
The manual_annotations table has to follow the specific layout given in the data folder in order to work with this module.
"""

__author__ = "Carolin Brune, Tobias Fehrenbach, Famke Baeuerle and Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################

from Bio.KEGG import REST
import cobra
import click
import yaml
import os
import re
import gffutils
import sqlalchemy
import logging
import pandas as pd

from ols_client import EBIClient
from Bio import Entrez, SeqIO
from refinegems.utility.databases import PATH_TO_DB, initialise_database
from libsbml import Model as libModel
from libsbml import SBMLReader, writeSBMLToFile, SBMLValidator, SBMLDocument
from datetime import date
from typing import Literal

################################################################################
# functions
################################################################################

# models
# ------

def load_model(modelpath: str|list[str], package:Literal['cobra','libsbml']) -> cobra.Model|list[cobra.Model]|libModel|list[libModel]:
    """Load a model. 

    Args:
        modelpath (str | list[str]): Path to the model or list of paths to models (string format).
        package (Literal['cobra','libsbml']): Package to use to load the model.

    Returns:
        cobra.Model|list[cobra.Model]|libModel|list[libModel]: The loaded model(s).
    """

    def load_cobra_model(modelpath:str) -> cobra.Model:
        """Load a model using COBRApy.

        Args:
            modelpath (str): Path to the model.
                Can be a xml, json, yml or mat file.

        Raises:
            ValueError: Unknown file extension

        Returns:
            cobra.Model: The loaded model object.
        """
        extension = os.path.splitext(modelpath)[1].replace('.','')

        match extension:
            case 'xml':
                data = cobra.io.read_sbml_model(modelpath)
            case 'json':
                data = cobra.io.load_json_model(modelpath)
            case 'yml':
                data = cobra.io.load_yaml_model(modelpath)
            case 'mat':
                data = cobra.io.load_matlab_model(modelpath)
            case _:
                raise ValueError('Unknown file extension for model: ', extension)

        return data
    
    def load_libsbml_model(modelpath:str) -> libModel:
        """Load a model with libsbml.

        Args:
            modelpath (str): Path to the model. Should be xml.

        Returns:
            libModel: The loaded model object
        """

        reader = SBMLReader()
        read = reader.readSBMLFromFile(modelpath)  # read from file
        mod = read.getModel()

        return mod

    match modelpath:
        # read in multiple models
        case list():

            loaded_models = []
            for modelpath in modelpath:
                if package == 'cobra':
                    loaded_models.append(load_cobra_model(modelpath))
                elif package == 'libsbml':
                    loaded_models.append(load_libsbml_model(modelpath))
            return loaded_models
        
        # read in a single model
        case str():

                if package == 'cobra':
                    return load_cobra_model(modelpath)
                elif package == 'libsbml':
                    return load_libsbml_model(modelpath)


def load_document_libsbml(modelpath: str) -> SBMLDocument:
    """Loads model document using libSBML

    Args:
        - modelpath (str): Path to GEM

    Returns:
        SBMLDocument: Loaded document by libSBML
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    return read


def write_model_to_file(model:libModel|cobra.Model, filename:str):
    """Save a model into a file.

    Args:
        - model (libModel|cobra.Model): 
            The model to be saved
        - filename (str): 
            The filename to save the model to.

    Raises:
        - ValueError: Unknown file extension for model
        - TypeError: Unknown model type
    """

    # save cobra model
    if isinstance(model, cobra.core.model.Model):
        try:
            extension = os.path.splitext(filename)[1].replace('.','')
            match extension:
                case 'xml':
                    cobra.io.write_sbml_model(model, filename)
                case 'json':
                    cobra.io.save_json_model(model, filename)
                case 'yml':
                    cobra.io.save_yaml_model(model, filename)
                case 'mat':
                    cobra.io.save_matlab_model(model, filename)
                case _:
                    raise ValueError('Unknown file extension for model: ', extension)
            logging.info("Modified model written to " + filename)
        except (OSError) as e:
            print("Could not write to file. Wrong path?")

    # save libsbml model
    elif isinstance(model, libModel):
        try:
            new_document = model.getSBMLDocument()
            writeSBMLToFile(new_document, filename)
            logging.info("Modified model written to " + filename)
        except (OSError) as e:
            print("Could not write to file. Wrong path?")
    # unknown model type or no model        
    else:
        message = f'Unknown model type {type(model)}. Cannot save.'
        raise TypeError(message)
    

# media
# -----

# other
# -----
# @TODO: sort more, make it more readable

def load_manual_annotations(tablepath: str='data/manual_curation.xlsx', sheet_name: str='metab') -> pd.DataFrame:
    """Loads metabolite sheet from manual curation table

    Args:
        - tablepath (str): Path to manual curation table. Defaults to 'data/manual_curation.xlsx'.
        - sheet_name (str): Sheet name for metabolite annotations. Defaults to 'metab'.

    Returns:
        pd.DataFrame: Table containing specified sheet from Excel file
    """
    man_ann = pd.read_excel(tablepath, sheet_name)
    return man_ann


def load_a_table_from_database(table_name_or_query: str, query: bool=True) -> pd.DataFrame:
    """| Loads the table for which the name is provided or a table containing all rows for which the query evaluates to 
       | true from the refineGEMs database ('data/database/data.db')

    Args:
        - table_name_or_query (str): Name of a table contained in the database 'data.db'/ a SQL query
        - query (bool): Specifies if a query or a table name is provided with table_name_or_query

    Returns:
        pd.DataFrame: Containing the table for which the name was provided from the database 'data.db'
    """
    table_name_or_query = sqlalchemy.text(table_name_or_query) if query else table_name_or_query
    sqlalchemy_engine_input = f'sqlite:///{PATH_TO_DB}'
    engine = sqlalchemy.create_engine(sqlalchemy_engine_input)
    open_con = engine.connect()

    db_table = pd.read_sql(table_name_or_query, open_con)

    open_con.close()
    return db_table


def load_manual_gapfill(tablepath: str='data/manual_curation.xlsx' , sheet_name: str='gapfill') -> pd.DataFrame:
    """Loads gapfill sheet from manual curation table

    Args:
        - tablepath (str): Path to manual curation table. Defaults to 'data/manual_curation.xlsx'.
        - sheet_name (str): Sheet name for reaction gapfilling. Defaults to 'gapfill'.

    Returns:
        pd.DataFrame: Table from Excel file sheet with name 'gapfill'/ specified sheet_name
    """
    man_gapf = pd.read_excel(tablepath, sheet_name)
    return man_gapf


def parse_dict_to_dataframe(str2list: dict) -> pd.DataFrame:
    """| Parses dictionary of form {str: list} & 
       | Transforms it into a table with a column containing the strings and a column containing the lists

    Args:
        str2list (dict): Dictionary mapping strings to lists

    Returns:
        pd.DataFrame: Table with column containing the strings and column containing the lists
    """
    # Get max number of list length
    max_len_of_list = max(map(len, str2list.values()))

    # Fill lists with None until all lists have the same size -> Required for pd.DataFrame
    for key in str2list:
       current_list = str2list.get(key)
       while len(current_list) != max_len_of_list:
          str2list.get(key).append(None)

    df = pd.DataFrame.from_dict(str2list).stack().T.reset_index()
    df = df.drop('level_0', axis=1)
    
    return df


def write_df_to_xlsx(dataframe: pd.DataFrame, filepath: str):
    """Writes reports stored in dataframes to xlsx file

    Args:
        - dataframe (pd.DataFrame): Table containing output
        - filepath (str): Path to file with filename
    """
    writer = pd.ExcelWriter(str(os.path.abspath('.')) + '/' + filepath)
    dataframe.to_excel(writer)
    writer.save()


def validate_libsbml_model(model: libModel) -> int:
    """Debug method: Validates a libSBML model with the libSBML validator
    
    Args:
        - model (libModel): A libSBML model
        
    Returns:
        int: Integer specifying if validate was successful or not
    """
    validator = SBMLValidator()
    doc = model.getSBMLDocument()
    
    return validator.validate(doc)


def parse_fasta_headers(filepath: str, id_for_model: bool=False) -> pd.DataFrame:
    """Parses FASTA file headers to obtain:
    
        - the protein_id
        - and the model_id (like it is obtained from CarveMe)
            
    corresponding to the locus_tag
        
    Args:
        - filepath (str): Path to FASTA file
        - id_for_model (bool): True if model_id similar to autogenerated GeneProduct ID should be contained in resulting table
        
    Returns:
        pd.DataFrame: Table containing the columns locus_tag, Protein_id & Model_id
    """
    keyword_list = ['protein', 'locus_tag']
    tmp_dict = dict()
    if id_for_model:
        locus2ids = {
            'locus_tag': [],
            'protein_id': [],
            'model_id': [],
            'name': []
        }
    else:
        locus2ids = {
            'locus_tag': [],
            'protein_id': [],
            'name': []
        }
   
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.description
            protein_id = record.id.split('|')[1].split('prot_')[1].split('.')[0].strip()
            descriptors = re.findall('\[+(.*?)\]', header)
            if id_for_model:
                model_id = re.sub("\||\.", "_", record.id)
                model_id = f'G_{model_id}'
         
            descriptors.insert(0, protein_id)
            
            tmp_dict['protein_id'] = str(protein_id)
            
            for entry in descriptors:
                entry = entry.strip().split('=')
               
                if entry[0] in keyword_list:
                    if entry[0] == 'protein_id':
                        tmp_dict[entry[0]] = entry[1].split('.')[0]
                    else:
                        tmp_dict[entry[0]] = entry[1]
                    
            locus2ids.get('locus_tag').append(tmp_dict.get('locus_tag'))
            locus2ids.get('protein_id').append(tmp_dict.get('protein_id'))
            locus2ids.get('name').append(tmp_dict.get('protein'))
            if id_for_model:
                locus2ids.get('model_id').append(model_id)
            
    return pd.DataFrame(locus2ids)


def search_ncbi_for_gpr(locus: str) -> str:
    """Fetches protein name from NCBI

    Args:
        - locus (str): NCBI compatible locus_tag

    Returns:
        str: Protein name|description
    """
    handle = Entrez.efetch(
        db="protein",
        id=locus,
        rettype="gbwithparts",
        retmode='text')
    records = SeqIO.parse(handle, "gb")

    for i, record in enumerate(records):
        if (locus[0] == 'W'):
            return record.description, locus
        else:
            for feature in record.features:
                if feature.type == "CDS":
                    return record.description, feature.qualifiers["locus_tag"][0]


def parse_gff_for_refseq_info(gff_file: str) -> pd.DataFrame:
    """Parses the RefSeq GFF file to obtain a mapping from the locus tag to the corresponding RefSeq identifier

    Args:
        gff_file (str): RefSeq GFF file of the input organism

    Returns:
        pd.DataFrame: Table mapping locus tags to their respective RefSeq identifiers
    """

    locus_tag2id = {}
    locus_tag2id['LocusTag'] = []
    locus_tag2id['ProteinID'] = []
      
    gff_db = gffutils.create_db(gff_file, ':memory:', merge_strategy='create_unique')

    for feature in gff_db.all_features():
        
        if (feature.featuretype == 'gene') and ('old_locus_tag' in feature.attributes):  # Get locus_tag & old_locus_tag
            current_locus_tag = feature.attributes['locus_tag']
            locus_tag2id['LocusTag'].append(feature.attributes['old_locus_tag'][0])
        elif (feature.featuretype == 'gene') and ('locus_tag' in feature.attributes):
            current_locus_tag = feature.attributes['locus_tag']
            locus_tag2id['LocusTag'].append(feature.attributes['locus_tag'][0])
            
        if (feature.featuretype == 'CDS') and ('protein_id' in feature.attributes): # Check if CDS has protein_id
            if feature.attributes['locus_tag'] == current_locus_tag:# Get protein_id if locus_tag the same
                locus_tag2id['ProteinID'].append(feature.attributes['protein_id'][0])

    locus_tag2id['LocusTag'] = locus_tag2id.get('LocusTag')[:len(locus_tag2id.get('ProteinID'))]

    return pd.DataFrame(locus_tag2id)


def parse_gff_for_gp_info(gff_file: str) -> pd.DataFrame:
    """Parses gff file of organism to find gene protein reactions based on locus tags

    Args:
        - gff_file (str): Path to gff file of organism of interest

    Returns:
        pd.DataFrame: Table containing mapping from locus tag to GPR
    """
    db = gffutils.create_db(
        gff_file,
        ':memory:',
        merge_strategy='create_unique')
    mapping_cds = {}
    for feature in db.all_features():
        attr = dict(feature.attributes)
        try:
            if str(attr['gbkey'][0]) == 'CDS':
                mapping_cds[attr['Name'][0]] = attr['Parent'][0]
        except BaseException:
            pass
    mapping_df = pd.DataFrame.from_dict(
        mapping_cds,
        columns=['Parent'],
        orient='index').reset_index().rename(
        columns={
            'index': 'GPR'})

    def extract_locus(feature):
        try:
            return db[feature].attributes['old_locus_tag'][0]
        except BaseException:
            pass
        return None

    mapping_df['locus_tag'] = mapping_df.apply(
        lambda row: extract_locus(row['Parent']), axis=1)
    return mapping_df.drop('Parent', axis=1)


def search_sbo_label(sbo_number: str) -> str:
    """Looks up the SBO label corresponding to a given SBO Term number

    Args:
        - sbo_number (str): Last three digits of SBO-Term as str

    Returns:
        str: Denoted label for given SBO Term
    """
    sbo_number = str(sbo_number)
    client = EBIClient()
    sbo = client.get_term('sbo', 'http://biomodels.net/SBO/SBO_0000' + sbo_number)
    return sbo['_embedded']['terms'][0]['label']


def save_user_input(configpath: str) -> dict[str: str]:
    """This aims to collect user input from the command line to create a config file, 
    will also save the user input to a config if no config was given

    Args:
        - configpath (str): Path to config file if present
        
    Returns:
        dict: Either loaded config file or created from user input
    """
    if os.path.isfile(configpath):
        with open(configpath) as f:
            config = yaml.safe_load(f)
        print(config)
        return config
    else:
        print('No config or no valid config given, you will be asked for input')
        user_input = {}
        
        update_db = click.confirm('Do you want to update the database?')
        user_input['db_update'] = update_db
        
        if update_db:
            initialise_database()

        out_path = click.confirm('Do you want to keep the output path "../rg_out/"?', default=True)
        if not out_path:
            user_input['out_path'] = click.prompt('Enter your desired output path')
        else:
            user_input['out_path'] = '../rg_out/'
        
        user_input['visualize'] = click.confirm('Do you want to generate visualizations of your model(s)?')
            
        growth_basis = click.prompt('Enter the base uptakes for growth simulation (d for default_uptake, m for minimal_uptake)')
        if growth_basis == 'd':
            user_input['growth_basis'] = 'default_uptake'
        if growth_basis == 'm':
            user_input['growth_basis'] = 'minimal_uptake'
            
        user_input['anaerobic_growth'] = click.confirm('Do you want to simulate anaerobic growth?')
        
        multiple = click.confirm('Do you want to simulate and compare multiple models?')
        user_input['multiple'] = multiple
        if multiple:
            list_of_models = []
            while True:
                file_path = click.prompt('Enter file path to model (or "stop" to stop)')
                if file_path.lower() == 'stop':
                    break
                elif os.path.isfile(file_path):
                    list_of_models.append(file_path)
                    print('Added file:', file_path)
                else:
                    print('File does not exist. Please enter a valid file path.')
            print('The following models will be compared:')
            print(list_of_models)
            user_input['multiple_paths'] = list_of_models
        possible_media = load_a_table_from_database('medium', False)['name'].to_list()
        possible_media_str = '|'.join(possible_media)
        list_of_media = []
        while True:
            medium = click.prompt(f'Enter medium to simulate growth on ({possible_media_str}) (or "stop" to stop)')
            if medium.lower() == 'stop':
                break
            elif medium in possible_media:
                if medium not in list_of_media:
                    list_of_media.append(medium)
                else:
                    print(medium + ' is already in the list.')
            else:
                print('Please choose a medium from the given list.')
        user_input['media'] = list_of_media
        
        single = click.confirm('Do you want to investigate or curate a single model?')
        user_input['single'] = single
        if single:
            not_valid = True
            while not_valid:
                model = click.prompt('Path to your model file.')
                if os.path.isfile(model):
                    user_input['model'] = model
                    not_valid = False
                else:
                    print('File does not exist. Please enter a valid file path')
            user_input['memote'] = click.confirm('Do you want to run MEMOTE (takes some time)?')    
            user_input['modelseed'] = click.confirm('Do you want to compare your model entities to the ModelSEED database?')
        
            gap_analysis = click.confirm('Do you want to run a gap analysis?') 
            user_input['gap_analysis'] = gap_analysis
            if gap_analysis:
                gap_analysis_params = {}
                db_to_compare = click.prompt('One of the choices KEGG|BioCyc|KEGG+BioCyc') #|GFF
                gap_analysis_params['db_to_compare'] = db_to_compare
                if db_to_compare == 'KEGG' or db_to_compare == 'KEGG+BioCyc':
                    user_input['organismid'] = click.prompt('Enter the KEGG Organism code')
                    user_input['gff_file'] = click.prompt('Enter the path to your organisms RefSeq GFF file')
                if db_to_compare == 'BioCyc' or db_to_compare == 'KEGG+BioCyc':
                    Path0 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with the columns \'Accession-2\' and \'Reaction of gene\'')
                    Path1 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with all reaction relevant information')
                    Path2 = click.prompt('Enter the path to your Biocyc TXT file containing a SmartTable with all metabolite relevant information')
                    Path3 = click.prompt('Enter path to protein FASTA file used as input for CarveMe')
                    gap_analysis_params['biocyc_files'] = [Path0, Path1, Path2, Path3]
                user_input['gap_analysis_params'] = gap_analysis_params
                
            mod = click.confirm('Do you want to use functions to modify your model?')
            if mod:
                
                new_path = click.confirm('Do you want to save your modified model to ' + user_input['out_path'] + '<model.id>_modified_<today>.xml?')
                if new_path:
                    user_input['model_out'] = 'stdout'
                else:
                    user_input['model_out'] = click.prompt('Enter path and filename to where to save the modified model')
                
                gapfill_model = click.confirm('Do you want to fill gaps in your model?')
                user_input['gapfill_model'] = gapfill_model
                
                if gapfill_model:
                    if not gap_analysis:
                        user_input['gap_analysis_file'] = click.prompt('Enter path to Excel file with which gaps should be filled')
                
                user_input['keggpathways'] = click.confirm('Do you want to add KEGG Pathways?')
                    
                user_input['sboterms'] = click.confirm('Do you want to update the SBO Terms?')
                
                user_input['charge_corr'] = click.confirm('Do you want to add charges to uncharged metabolites?')
                    
                man_cur = click.confirm('Do you want to modify your model with the manual curations table?')
                user_input['man_cur'] = man_cur

                if man_cur:
                    entrez_email = click.prompt('Email to access NCBI Entrez')
                    user_input['entrez_email'] = entrez_email
                    man_cur_type = click.prompt('Enter type of curation (gapfill|metabs)')
                    user_input['man_cur_type'] = man_cur_type
                    man_cur_table = click.prompt('Enter the path to the manual curations table')
                    user_input['man_cur_table'] = man_cur_table

                polish = click.confirm('Do you want to polish the model?')
                user_input['polish'] = polish

                if polish:
                    entrez_email = click.prompt('Email to access NCBI Entrez')
                    user_input['entrez_email'] = entrez_email
                    id_db = click.prompt('What database is your model based on? BIGG|VMH')
                    user_input['id_db'] = id_db
                    gff_file = click.prompt('If possible, provide the path to the RefSeq GFF file of your organism')
                    user_input['gff_file'] = gff_file if gff_file != 'None' else None
                    protein_fasta = click.prompt('If possible, provide the path to the Protein FASTA file used for CarveMe')
                    user_input['protein_fasta'] = protein_fasta if protein_fasta != 'None' else None
                    lab_strain = not click.confirm('Does your modeled organism have a database entry?', default=True)
                    user_input['lab_strain'] = lab_strain
                    organismid = click.prompt('If possible, provide the KEGG organism code of your organism')
                    user_input['organismid'] = organismid if (organismid != 'None') else None
                    
                biomass = click.confirm('Do you want to check & normalise the biomass function(s)?')
                user_input['biomass'] = biomass
                    
            else:
                user_input['keggpathways'] = False
                user_input['polish'] = False
                user_input['biomass'] = False
                user_input['sboterms'] = False
                user_input['charge_corr'] = False
                user_input['gapfill_model'] = False
                user_input['man_cur'] = False
            
        today = date.today().strftime("%Y%m%d")
        
        print('This is your input:')
        print(user_input)
        if not os.path.isdir(user_input['out_path']):
            print('Given out_path is not yet a directory, creating ' + user_input['out_path'])
            os.makedirs(user_input['out_path'])
        with open(user_input['out_path'] + 'user_input_' + str(today) + '.yaml', 'w') as f:
            yaml.dump(user_input, f)
        print('Your input was saved as yaml to '+ user_input['out_path'] + 'user_input_' + str(today) + '.yaml')
        return user_input


# KEGG
# ----
# @TODO
def kegg_reaction_parser(rn_id):
    """Get the KEGG reaction entry for a KEGG reaction ID
    and parse it to retrieve information in form of a dictionary.

    :param rn_id: The KEGG reaction ID that shall be retrived and parsed.
    :type  rn_id: string
    :returns:     Information about name, equationand database links of the entry (reaction).
    :rtype:       dict, keys include 'name','db','rc' and 'equation'
    """

    # get KEGG reaction entry
    # .................
    #@TODO
    #     implement try-catch for errors (just in case)
    #     although they (theoreticlly) should not occur
    # current solution: do it outside this function
    # .................
    kegg_reac = REST.kegg_get(F'rn:{rn_id}')
    kegg_reac = kegg_reac.read()

    # parse the entry for necessary information
    features = {}
    db_entries = []
    pathways = []
    rc = []
    references = {}
    collect = False
    for line in kegg_reac.split('\n'):
        if line:
            if line.startswith('NAME'):
                features['name'] = line.replace('NAME','',1).strip()
            elif line.startswith('EQUATION'):
                features['equation'] = line.replace('EQUATION','',1).strip()
            elif line.startswith('ENZYME'):
                references['ec-code'] = line.replace('ENZYME','',1).strip()
            elif line.startswith('RCLASS'):
                rc.append(line.replace('RCLASS','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('PATHWAY'):
                pathways.append(line.replace('PATHWAY','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('DBLINKS'):
                db_entries.append(line.replace('DBLINKS','',1).strip())
                collect = True
            elif collect == True and line[0] != '/':
                if len(db_entries) == 0:
                    if line[0].isupper():
                        collect = False
                    else:
                        line = line.strip()
                        if line.startswith('RC'):
                            rc.append(line.split(' ')[0])
                        else:
                            pathways.append(line.split(' ')[0])
                else:
                    db_entries.append(line.strip())
            else:
                continue

    # parse references
    for entry in db_entries:
        db, identifier = entry.split(':')
        db = db.strip().lower()
        if db in references:
            references[db] = references[db].append(identifier)
        else:
            references[db] = [identifier.strip()]
    if len(pathways) > 0:
        references['kegg.pathway'] = pathways

    if len(references) > 0:
        features['db'] = references
    if len(rc) > 0:
        features['rc'] = rc

    return features
