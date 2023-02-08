#!/usr/bin/env python
import re
import gffutils
import pandas as pd
from Bio import Entrez, SeqIO

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"


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
        'Protein_id': [],
        'Model_id': []
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
            locus2ids.get('Protein_id').append(tmp_dict.get('protein_id'))
            locus2ids.get('Model_id').append(model_id)
            
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
 