#!/usr/bin/env python
""" Provides functions to load and write models, media definitions and the manual annotation table

Depending on the application the model needs to be loaded with cobra (memote)
or with libSBML (activation of groups). The media definitions are denoted in a csv within the data folder of this repository, thus the functions will only work if the user clones the repository. The manual_annotations table has to follow the specific layout given in the data folder in order to work with this module.
"""

import cobra
import click
import yaml
import os
import re
import gffutils
import sqlalchemy
import pandas as pd
from ols_client import EBIClient
from Bio import Entrez, SeqIO
from refinegems.databases import PATH_TO_DB
from libsbml import SBMLReader, writeSBMLToFile, Model, SBMLValidator, SBMLDocument
from datetime import date

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


def load_model_cobra(modelpath: str) -> cobra.Model:
    """loads model using cobrapy

    Args:
        modelpath (str): Path to GEM

    Returns:
        cobra-model: loaded model by cobrapy
    """
    mod = cobra.io.read_sbml_model(modelpath)
    return mod


def load_model_libsbml(modelpath: str) -> Model:
    """loads model using libsbml

    Args:
        modelpath (str): Path to GEM

    Returns:
        libsbml-model: loaded model by libsbml
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    mod = read.getModel()
    return mod


def load_multiple_models(models: list[str], package: str) -> list:
    loaded_models = []
    for modelpath in models:
        if package == 'cobra':
            loaded_models.append(load_model_cobra(modelpath))
        elif package == 'libsbml':
            loaded_models.append(load_model_libsbml(modelpath))
    return loaded_models


def load_document_libsbml(modelpath: str) -> SBMLDocument:
    """loads model document using libsbml

    Args:
        modelpath (str): Path to GEM

    Returns:
        libsbml-document: loaded document by libsbml
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    return read


def load_medium_custom(mediumpath: str) -> pd.DataFrame:
    """Helper function to read medium csv

    Args:
        mediumpath (Str): path to csv file with medium

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium['BiGG_R'] = 'R_EX_' + medium['BiGG'] + '_e'
    medium['BiGG_EX'] = 'EX_' + medium['BiGG'] + '_e'
    return medium


def load_medium_from_db(mediumname: str) -> pd.DataFrame:
    """Wrapper function to extract subtable for the requested medium from the database 'data.db'

    Args:
        mediumname (Str): name of medium to test growth on

    Returns:
        df: pandas dataframe containing composition for one medium with metabs added as BiGG_EX exchange reactions
    """
    medium_query = f"SELECT * FROM media m JOIN media_compositions mc ON m.id = mc.medium_id WHERE m.medium = '{mediumname}'" 
    medium = load_a_table_from_database(medium_query)
    medium = medium[['medium', 'medium_description', 'BiGG', 'substance']]
    medium['BiGG_R'] = 'R_EX_' + medium['BiGG'] + '_e'
    medium['BiGG_EX'] = 'EX_' + medium['BiGG'] + '_e'
    return medium


def load_all_media_from_db(mediumpath: str) -> pd.DataFrame: 
    """Helper function to extract media definitions from media_db.csv

    Args:
        mediumpath (Str): path to csv file with medium database

    Returns:
        df: pandas dataframe of csv with metabs added as BiGG_EX exchange reactions
    """
    media = pd.read_csv(mediumpath, sep=';')
    media['BiGG_R'] = 'R_EX_' + media['BiGG'] + '_e'
    media['BiGG_EX'] = 'EX_' + media['BiGG'] + '_e'

    media['group'] = media['medium'].ne(media['medium'].shift()).cumsum()
    grouped = media.groupby('group')
    media_dfs = []
    for name, data in grouped:
        media_dfs.append(data.reset_index(drop=True))
    return media_dfs


def load_manual_annotations(tablepath: str='data/manual_curation.xlsx', sheet_name: str='metab') -> pd.DataFrame:
    """loads metabolite sheet from manual curation table

    Args:
        tablepath (str): Path to manual curation table. Defaults to 'data/manual_curation.xlsx'.
        sheet_name (str): Sheet name for metabolite annotations. Defaults to 'metab'.

    Returns:
        df: table as pandas df
    """
    man_ann = pd.read_excel(tablepath, sheet_name)
    return man_ann


def load_a_table_from_database(table_name_or_query: str) -> pd.DataFrame:
    """Loads the table for which the name is provided or a table containing all rows for which the query evaluates to 
        true from the refineGEMs database ('data/database/data.db')

    Args:
        table_name_or_query (str): Name of a table contained in the database 'data.db'/ a SQL query

    Returns:
        pd.DataFrame: Containing the table for which the name was provided from the database 'data.db'
    """
    sqlalchemy_engine_input = f'sqlite:///{PATH_TO_DB}'
    engine = sqlalchemy.create_engine(sqlalchemy_engine_input)
    open_con = engine.connect()
    
    db_table = pd.read_sql(table_name_or_query, open_con)
    
    open_con.close()
    return db_table


def load_manual_gapfill(tablepath: str='data/manual_curation.xlsx' , sheet_name: str='gapfill') -> pd.DataFrame:
    """loads gapfill sheet from manual curation table

    Args:
        tablepath (str): Path to manual curation table. Defaults to 'data/manual_curation.xlsx'.
        sheet_name (str): Sheet name for reaction gapfilling. Defaults to 'gapfill'.

    Returns:
        df: table as pandas df
    """
    man_gapf = pd.read_excel(tablepath, sheet_name)
    return man_gapf


def write_to_file(model: Model, new_filename: str):
    """Writes modified model to new file

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename / path for modified model
    """
    new_document = model.getSBMLDocument()
    writeSBMLToFile(new_document, new_filename)
    print("Modified model written to " + new_filename)


def write_report(dataframe: pd.DataFrame, filepath: str):
    """Writes reports stored in dataframes to xlsx file

    Args:
        dataframe (pd.DataFrame): table containing output
        filepath (str): path to file with filename
    """
    writer = pd.ExcelWriter(str(os.path.abspath('.')) + '/' + filepath)
    dataframe.to_excel(writer)
    writer.save()


def validate_libsbml_model(model: Model):
    ''' Debug method: Validates a libSBML model with the libSBML validator
    
        Params:
            - model (Model): A libSBML Model
    '''
    validator = SBMLValidator()
    doc = model.getSBMLDocument()
    
    return validator.validate(doc)


def parse_fasta_headers(filepath: str, id_for_model: bool=False) -> pd.DataFrame:
    """Parses FASTA file headers to obtain:
        - the protein_id
        - and the model_id (like it is obtained from CarveMe)
        corresponding to the locus_tag
        
        Args:
            filepath (str): Path to FASTA file
            
        Returns:
            a pandas dataframe containing the columns locus_tag, Protein_id & Model_id
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


def search_ncbi_for_gpr(locus):
    """fetches protein name from NCBI

    Args:
        locus (string): NCBI compatible locus_tag

    Returns:
        str: protein name / description
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


def parse_gff_for_gp_info(gff_file):
    """Parses gff file of organism to find gene protein reactions based on locus tags

    Args:
        gff_file (Str): path to gff file of organism of interest

    Returns:
        df: table containing mapping from locus tag to GPR
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
    """looks up the SBO label corresponding to a given SBO Term number

    Args:
        sbo_number (str): Last three digits of SBO-Term as str

    Returns:
        str: denoted label for given SBO Term
    """
    sbo_number = str(sbo_number)
    client = EBIClient()
    sbo = client.get_term('sbo', 'http://biomodels.net/SBO/SBO_0000' + sbo_number)
    return sbo['_embedded']['terms'][0]['label']


def save_user_input(configpath):
    """This aims to collect user input from the command line to create a config file, 
    will also save the user input to a config if no config was given

    Args:
        configpath (str): path to config file if present
        
    Returns:
        dict: either loaded config file or created from user input
    """
    if os.path.isfile(configpath):
        with open(configpath) as f:
            config = yaml.safe_load(f)
        print(config)
        return config
    else:
        print('No config or no valid config given, you will be asked for input')
        user_input = {}
        not_valid = True
        while not_valid:
            model = click.prompt('Path to your model file.')
            if os.path.isfile(model):
                user_input['model'] = model
                not_valid = False
            else:
                print('File does not exist. Please enter a valid file path')

        out_path = click.confirm('Do you want to keep the out path "../rg_out/"?', default=True)
        if not out_path:
            user_input['out_path'] = click.prompt('Enter you desired output path')
        else:
            user_input['out_path'] = '../rg_out/'
            
        growth_basis = click.prompt('Enter the base uptakes for growth simulation (d for default_uptake, m for minimal_uptake)')
        if growth_basis == 'd':
            user_input['growth_basis'] = 'default_uptake'
        if growth_basis == 'm':
            user_input['growth_basis'] = 'minimal_uptake'
        
        multiple = click.confirm('Do you want to simulate and compare multiple models?')
        user_input['multiple'] = multiple
        if multiple:
            user_input['visualize'] = True
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
        list_of_media = []
        while True:
            medium = click.prompt('Enter medium to simulate growth on (SNM3|LB|M9|SMM|CGXII|RPMI) (or "stop" to stop)')
            if medium.lower() == 'stop':
                break
            elif medium in ['SNM3', 'RPMI', 'CGXlab', 'LB', 'M9', 'CGXII', 'CasA']:
                if medium not in list_of_media:
                    list_of_media.append(medium)
                else:
                    print(medium + ' is already in the list.')
            else:
                print('Please choose a medium from the given list.')
        user_input['media'] = list_of_media
        user_input['memote'] = click.confirm('Do you want to run MEMOTE (takes some time)?')    
        user_input['modelseed'] = click.confirm('Do you want to compare your model entities to the ModelSEED database?')
        user_input['output'] = 'xlsx'
        
        gapfill_analysis = click.confirm('Do you want to run the gapfill analysis?') 
        user_input['gapfill_analysis'] = gapfill_analysis
        if gapfill_analysis:
            gapfill_params = {}
            db_to_compare = click.prompt('One of the choices KEGG|BioCyc|GFF|KEGG+BioCyc')
            gapfill_params['db_to_compare'] = db_to_compare
            if db_to_compare == 'KEGG' or db_to_compare == 'KEGG+BioCyc':
                gapfill_params['organismid'] = click.prompt('Enter the KEGG Organism ID')
            if db_to_compare == 'GFF':
                gapfill_params['gff_file'] = click.prompt('Enter the path to your organisms GFF file')
            if db_to_compare == 'BioCyc' or db_to_compare == 'KEGG+BioCyc':
                Path0 = click.prompt()
                Path1 = click.prompt()
                Path2 = click.prompt()
                Path3 = click.prompt()
                gapfill_params['biocyc_files'] = [Path0, Path1, Path2, Path3]
            user_input['gapfill_analysis_params'] = gapfill_params
        else:
            user_input['gapfill_model'] = False
            
        mod = click.confirm('Do you want to use functions to modify your model?')
        if mod:
            if gapfill_analysis:
                user_input['gapfill_model'] = click.confirm('Do you want to gap fill your model?')
            
            kegg = click.confirm('Do you want to add KEGG Pathways?')
            user_input['keggpathways'] = kegg

            if kegg:
                kegg_path = click.prompt('Enter the modified file name')
                user_input['kegg_path'] = kegg_path
                
            polish = click.confirm('Do you want to polish the model?')
            user_input['polish'] = polish

            if polish:
                entrez_email = click.prompt('Email to access NCBI Entrez')
                user_input['entrez_email'] = entrez_email
                id_db = click.prompt('What database is your model based on? BIGG|VMH')
                user_input['id_db'] = id_db
                lab_strain = click.confirm('Does your modeled organism have NO database entry?')
                user_input['lab_strain'] = lab_strain
                protein_fasta = click.prompt('If possible, provide the path to your Protein FASTA file used for CarveMe')
                user_input['protein_fasta'] = protein_fasta
                polish_path = click.prompt('Enter the modified file name')
                user_input['polish_path'] = polish_path
                
            sboterms = click.confirm('Do you want to update the SBO Terms?')
            user_input['sboterms'] = sboterms

            if sboterms:
                sbo_path = click.prompt('Enter the modified file name')
                user_input['sbo_path'] = sbo_path
            
            charge_corr = click.confirm('Do you want to add charges to uncharged metabolites?')
            user_input['charge_corr'] = charge_corr

            if charge_corr:
                charge_path = click.prompt('Enter the modified file name')
                user_input['charge_path'] = charge_path
                user_input['charge_report_path'] = '../rg_out/multiple_charges.csv'
                
            man_cur = click.confirm('Do you want to modify your model with the manual curations table?')
            user_input['man_cur'] = man_cur

            if man_cur:
                entrez_email = click.prompt('Email to access NCBI Entrez')
                user_input['entrez_email'] = entrez_email
                man_cur_type = click.prompt('Enter type of curation (gapfill|metabs)')
                user_input['man_cur_type'] = man_cur_type
                man_cur_table = click.prompt('Enter the path to the manual curations table')
                user_input['man_cur_table'] = man_cur_table
                man_cur_path = click.prompt('Enter the modified file name')
                user_input['man_cur_path'] = man_cur_path

        else:
            user_input['keggpathways'] = False
            user_input['polish'] = False
            user_input['sboterms'] = False
            user_input['charge_corr'] = False
            user_input['man_cur'] = False
            
        today = date.today().strftime("%Y%m%d")
        
        print('This is your input:')
        print(user_input)
        if not os.path.isdir(user_input['out_path']):
            print('Given out_path is not yet a directory, creating ' + user_input['out_path'])
            os.makedirs(user_input['out_path'])
        with open(user_input['out_path'] + 'user_input.yaml', 'w') as f:
            yaml.dump(user_input, f)
        print('Your input was saved as yaml to '+ user_input['out_path'] + 'user_input' + str(today) + '.yaml')
        return user_input