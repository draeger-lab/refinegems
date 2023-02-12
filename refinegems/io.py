#!/usr/bin/env python
""" Provides functions to load and write models, media definitions and the manual annotation table

Depending on the application the model needs to be loaded with cobra (memote)
or with libSBML (activation of groups). The media definitions are denoted in a csv within the data folder of this repository, thus the functions will only work if the user clones the repository. The manual_annotations table has to follow the specific layout given in the data folder in order to work with this module.
"""

import cobra
import os
import re
import gffutils
import pandas as pd
from Bio import Entrez, SeqIO
from libsbml import SBMLReader, writeSBMLToFile, Model, SBMLValidator, SBMLDocument

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


def load_medium_from_db(mediumpath: str, mediumname: str) -> pd.DataFrame:
    """Helper function to read standard media_db.csv

    Args:
        mediumpath (Str): path to csv file with medium database
        mediumname (Str): name of medium to test growth on

    Returns:
        df: pandas dataframe of csv
    """
    medium = pd.read_csv(mediumpath, sep=';')
    medium = medium.loc[medium['medium'] == mediumname]
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


def parse_fasta_headers_to_dict(filepath: str) -> dict[str: tuple[str, str]]:
    """Parses a .fasta file to reveal the protein IDs
       NCBI Protein headers contain information on locus tag and protein name
    
    Params:
        - filepath (str): path to .fasta file

    Returns:
        dict: mapping from protein ids to (protein name, locus tag)
    """
    keyword_list = ['protein', 'locus_tag']
    tmp_dict = dict()
    id2locus_names = dict()
    
    with open(filepath, 'r') as handle:
        
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.description
            identifier = record.id.split('|')[1].split('prot_')[1].split('.')[0].strip()
            descriptors = re.findall('\[+(.*?)\]',header)
            
            descriptors.insert(0, identifier)
                
            tmp_dict['protein_id'] = identifier
                
            for entry in descriptors:
                
                entry = entry.strip().split('=')
                
                if entry[0] in keyword_list:
                    if entry[0] == 'protein_id':
                        tmp_dict[entry[0]] = entry[1].split('.')[0]
                    else:
                        tmp_dict[entry[0]] = entry[1]
                
            id2locus_names[tmp_dict.get('protein_id')] = (tmp_dict.get('protein'), tmp_dict.get('locus_tag'))
                
    return id2locus_names

def parse_fasta_headers(filepath: str) -> pd.DataFrame:
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
    locus2ids = {
        'locus_tag': [],
        'protein_id': [],
        'model_id': [],
        'name': []
    }
   
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.description
            protein_id = record.id.split('|')[1].split('prot_')[1].split('.')[0]
            descriptors = re.findall('\[+(.*?)\]', header)
            model_id = re.sub("\||\.", "_", record.id)
            model_id = f'G_{model_id}'
         
            descriptors.insert(0, protein_id)
            
            tmp_dict['protein_id'] = protein_id
            
            for entry in descriptors:
                entry = entry.split('=')
               
                if entry[0] in keyword_list:
                    if entry[0] == 'protein_id':
                        tmp_dict[entry[0]] = entry[1].split('.')[0]
                    else:
                        tmp_dict[entry[0]] = entry[1]
                    
            locus2ids.get('locus_tag').append(tmp_dict.get('locus_tag'))
            locus2ids.get('protein_id').append(tmp_dict.get('protein_id'))
            locus2ids.get('name').append(tmp_dict.get('protein'))
            locus2ids.get('model_id').append(model_id)
            
    return pd.DataFrame(locus2ids)
 

# Function originally from analysis_kegg
def get_name_from_locus(locus):
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
        return record.description
 

# Function originally from analysis_kegg
def get_locus_gpr(gff_file):
    """Searches gff file of organism for gene protein reactions based on locus tags

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
 