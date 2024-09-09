#!/usr/bin/env python
""" Provides functions to load and write models, media definitions and the manual annotation table

Depending on the application the model needs to be loaded with COBRApy (e.g. memote)
or with libSBML (e.g. activation of groups). Some might even require both (e.g. gap filling).
The manual_annotations table has to follow the specific layout given in the data folder in order to work with this module.
"""

__author__ = "Carolin Brune, Tobias Fehrenbach, Famke Baeuerle and Gwendolyn O. Döbel"

################################################################################
# requirements
################################################################################


import cobra
import os
import re
import gffutils
import sqlalchemy
import logging
import pandas as pd

from ols_client import EBIClient
from Bio import SeqIO
from libsbml import Model as libModel
from libsbml import SBMLReader, writeSBMLToFile, SBMLValidator, SBMLDocument
from pathlib import Path
from typing import Literal, Union

from .databases import PATH_TO_DB

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# models
# ------

def load_model(modelpath: Union[str,list[str]], package:Literal['cobra','libsbml']) -> Union[cobra.Model,list[cobra.Model],libModel,list[libModel]]:
    """Load a model. 

    Args:
        - modelpath (str | list[str]): 
            Path to the model or list of paths to models (string format).
        - package (Literal['cobra','libsbml']): 
            Package to use to load the model.

    Returns:
        cobra.Model|list[cobra.Model]|libModel|list[libModel]: 
            The loaded model(s).
    """

    def load_cobra_model(modelpath:str) -> cobra.Model:
        """Load a model using COBRApy.

        Args:
            modelpath (str): Path to the model.
                Can be a xml, json, yml or mat file.

        Raises:
            - ValueError: Unknown file extension

        Returns:
            cobra.Model: 
                The loaded model object.
        """
        extension = os.path.splitext(modelpath)[1].replace('.','')

        match extension:
            case 'xml' | 'sbml':
                data = cobra.io.read_sbml_model(modelpath)
            case 'json':
                data = cobra.io.load_json_model(modelpath)
            case 'yml' | 'yaml':
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
            libModel: 
                The loaded model object
        """

        reader = SBMLReader()
        read = reader.readSBMLFromFile(modelpath)  # read from file
        mod = read.getModel()

        return mod

    match modelpath:
        # read in multiple models
        case list():

            loaded_models = []
            for m in modelpath:
                if package == 'cobra':
                    loaded_models.append(load_cobra_model(m))
                elif package == 'libsbml':
                    loaded_models.append(load_libsbml_model(m))
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
        - modelpath (str): 
            Path to GEM

    Returns:
        SBMLDocument: 
            Loaded document by libSBML
    """
    reader = SBMLReader()
    read = reader.readSBMLFromFile(modelpath)  # read from file
    return read


def write_model_to_file(model:Union[libModel,cobra.Model], filename:str):
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
    

# other
# -----


def load_a_table_from_database(table_name_or_query: str, query: bool=True) -> pd.DataFrame:
    """| Loads the table for which the name is provided or a table containing all rows for which the query evaluates to 
       | true from the refineGEMs database ('data/database/data.db')

    Args:
        - table_name_or_query (str): 
            Name of a table contained in the database 'data.db'/ a SQL query
        - query (bool): 
            Specifies if a query or a table name is provided with table_name_or_query

    Returns:
        pd.DataFrame: 
            Containing the table for which the name was provided from the database 'data.db'
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
        - tablepath (str): 
            Path to manual curation table. Defaults to 'data/manual_curation.xlsx'.
        - sheet_name (str): 
            Sheet name for reaction gapfilling. Defaults to 'gapfill'.

    Returns:
        pd.DataFrame: 
            Table from Excel file sheet with name 'gapfill'/ specified sheet_name
    """
    man_gapf = pd.read_excel(tablepath, sheet_name)
    return man_gapf


def parse_dict_to_dataframe(str2list: dict) -> pd.DataFrame:
    """| Parses dictionary of form {str: list} & 
       | Transforms it into a table with a column containing the strings and a column containing the lists

    Args:
        str2list (dict): 
            Dictionary mapping strings to lists

    Returns:
        pd.DataFrame: 
            Table with column containing the strings and column containing the lists
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


def validate_libsbml_model(model: libModel) -> int:
    """Debug method: Validates a libSBML model with the libSBML validator
    
    Args:
        - model (libModel): 
            A libSBML model
        
    Returns:
        int: 
            Integer specifying if validate was successful or not
    """
    validator = SBMLValidator()
    doc = model.getSBMLDocument()
    
    return validator.validate(doc)


# FASTA
# -----
# @DEPRECATE: This function could maybe be deprecated in a future update.
# @TODO: Check usage!
def parse_fasta_headers(filepath: str, id_for_model: bool=False) -> pd.DataFrame:
    """Parses FASTA file headers to obtain:
    
        - the protein_id
        - and the model_id (like it is obtained from CarveMe)
            
    corresponding to the locus_tag
        
    Args:
        - filepath (str): 
            Path to FASTA file
        - id_for_model (bool): 
            True if model_id similar to autogenerated GeneProduct ID should be contained in resulting table
        
    Returns:
        pd.DataFrame: 
            Table containing the columns locus_tag, Protein_id & Model_id
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


def create_missing_genes_protein_fasta(fasta,outdir, missing_genes):
    
    # format the missing genes' locus tags
    locus_tags_dict = {_:'[locus_tag='+_+']' for _ in list(missing_genes['locus_tag'])}
    
    # parse the protein FASTA
    protfasta = SeqIO.parse(fasta,'fasta')
    
    # extract the sequences of the missing genes only
    missing_seqs = []
    for seq in protfasta:
        # Case 1: locus tag equals FASTA header
        if seq.id in locus_tags_dict.keys():
            missing_seqs.append(seq)
        # Case 2: locus tag as descriptor
        elif any([k for k,v in locus_tags_dict.items() if v in seq.description]):
            for k,v in locus_tags_dict.items():
                if v in seq.description:
                    seq.id = k
                    missing_seqs.append(seq)
                    break
        # Case _: locus tag either in model or errounous
        else:
            pass
        
    # save the collected sequences in a new file
    if outdir:
        outfile = Path(outdir,'missing_genes.fasta')
    else:
        outfile = Path('missing_genes.fasta')
    SeqIO.write(missing_seqs, outfile, 'fasta')
    
    return outfile


# GFF
# ---

def parse_gff_for_refseq_info(gff_file: str) -> pd.DataFrame:
    """Parses the RefSeq GFF file to obtain a mapping from the locus tag to the corresponding RefSeq identifier

    Args:
        - gff_file (str): 
            RefSeq GFF file of the input organism

    Returns:
        - pd.DataFrame: 
            Table mapping locus tags to their respective RefSeq identifiers
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


def parse_gff_for_cds(gffpath, keep_attributes=None):
    # load the gff
    gff = gffutils.create_db(gffpath, ':memory:', merge_strategy="create_unique")
    # extract the attributes of the CDS 
    cds = pd.DataFrame.from_dict([_.attributes for _ in gff.features_of_type('CDS')])
    cds = cds.explode('locus_tag')
    genes = pd.DataFrame.from_dict([_.attributes for _ in gff.features_of_type('gene')])
    genes = genes.explode('locus_tag')
    if 'old_locus_tag' in genes.columns and 'locus_tag' in genes.columns:
        cds = cds.merge(genes[['locus_tag','old_locus_tag']], how='left', 
                        on='locus_tag')
        cds.drop(columns=['locus_tag'], inplace=True)
        cds.rename(columns={'old_locus_tag':'locus_tag'},inplace=True)
    # keep only certain columns
    if keep_attributes:
        cds = cds[[_ for _ in cds.columns if _ in list(keep_attributes.keys())]]
        # rename columns 
        cds.rename(columns={k:v for k,v in keep_attributes.items() if k in cds.columns}, inplace=True)
    if 'locus_tag' in cds.columns:
        cds = cds.explode('locus_tag')

    return cds


# else:
# -----

def search_sbo_label(sbo_number: str) -> str:
    """Looks up the SBO label corresponding to a given SBO Term number

    Args:
        - sbo_number (str): 
            Last three digits of SBO-Term as str

    Returns:
        str: 
            Denoted label for given SBO Term
    """
    sbo_number = str(sbo_number)
    client = EBIClient()
    sbo = client.get_term('sbo', 'http://biomodels.net/SBO/SBO_0000' + sbo_number)
    return sbo['_embedded']['terms'][0]['label']
