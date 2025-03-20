#!/usr/bin/env python
"""General functions to polish a model
"""

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel and Carolin Brune"
# @TODO Clean-up this module!
################################################################################
# requirements
################################################################################

import logging
import pandas as pd
import re

from Bio import Entrez

from datetime import date
from libsbml import GeneProduct

from tqdm.auto import tqdm

from ..utility.cvterms import add_cv_term_genes
from ..utility.db_access import search_ncbi_for_gpr
from ..utility.io import parse_fasta_headers

################################################################################
# functions
################################################################################

#--------------------------------- Function to add URIs from the IDs for GeneProducts ---------------------------------#
def cv_ncbiprotein(gene_list, email, locus2id: pd.DataFrame, protein_fasta: str, filename: str, lab_strain: bool=False):
    """Adds NCBI Id to genes as annotation

    Args:
        - gene_list (list): 
            libSBML ListOfGenes
        - email (str): 
            User Email to access the Entrez database
        - locus2id (pd.DataFrame): 
            Table mapping locus tags to their corresponding RefSeq/NCBI Protein identifiers
        - protein_fasta (str): 
            The path to the CarveMe protein.fasta input file
        - filename (str):
            Path to output file for genes where no annotation could be added
        - lab_strain (bool): 
            Needs to be set to True if strain was self-annotated
            and/or the locus tags in the CarveMe input file should be kept   
    """
    Entrez.email = email
    if locus2id is not None:
        locus2id = locus2id.groupby('ProteinID').agg(list)
                    
    id2locus_name = None  # Needs to be initialised, otherwise UnboundLocalError: local variable 'id2locus_name' referenced before assignment
    entry = None # Needs to be initialised, otherwise UnboundLocalError: local variable 'entry' referenced before assignment
    if protein_fasta:
        if protein_fasta.strip() != '': 
            id2locus_name = parse_fasta_headers(protein_fasta)
            id2locus_name.set_index('protein_id')
    
    genes_missing_annotation = []

    print('Setting CVTerms and removing notes for all genes:')
    for gene in tqdm(gene_list):
        # Get current gene_id
        gene_id = gene.getId()
        
        if not gene.isSetMetaId():
            gene.setMetaId('meta_' + gene_id)

        # Remove 'G_' prefix as for further handling unnecessary
        gene_id = gene_id.removeprefix('G_')
        
        if (gene_id[0] == 'W'): #addition to work with KC-Na-01
            # @TODO Implement better way of handling this -> Handle with COBRApy
            entry = gene_id[:-7] if '__46__' in gene_id else gene_id[:-2] # Required for VMH models
            add_cv_term_genes(entry, 'REFSEQ', gene)
            add_cv_term_genes(entry, 'NCBI', gene, lab_strain)
            name, locus = search_ncbi_for_gpr(entry)
            gene.setName(name)
            if (locus2id is not None) and (entry in locus2id.index):
                locus = locus2id.loc[entry, 'LocusTag']
            if len(locus) == 1: gene.setLabel(locus[0])
            else: genes_missing_annotation.append(gene_id)
        
        elif (gene_id != 'spontaneous') and (gene_id != 'Unknown'): # Has to be omitted as no additional data can be retrieved neither from NCBI nor the CarveMe input file
            if 'prot_' in gene_id:
                id_string = gene_id.split(r'prot_')[1].split(r'_')  # All NCBI CDS protein FASTA files have the NCBI protein identifier after 'prot_' in the FASTA identifier
                ncbi_id = id_string[0]  # If identifier contains no '_', this is full identifier
            else:
                id_string = gene_id.split(r'_')
                if 'peg' in id_string: 
                    genes_missing_annotation.append(gene_id)
                    continue
              
            if len(id_string) == 2: # Can be the case if ID is locus tag, for example
                genes_missing_annotation.append(gene_id)
                continue # Ignore locus tags as no valid identifiers
            if (len(id_string) > 2):  # Identifier contains '_'
            # Check that the second entry consists of a sequence of numbers -> Valid RefSeq identifier! 
            # (Needs to be changed if there are other gene idenitfiers used that could contain '_' & need to be handled differently)
                if re.fullmatch(r'^\d+\d+$', id_string[1], re.IGNORECASE):
                    # Merge the first two parts with '_' as this is complete identifier
                    # Merge the resulting string with the third string in i_string to get complete identifier with version spec
                    ncbi_id = f'{"_".join(id_string[:2])}.{id_string[2]}'
                else: # Can be the case if locus tag is part of the ID
                    genes_missing_annotation.append(gene_id)
                    continue # Ignore locus tags as no valid identifiers

            # If identifier matches RefSeq ID pattern   
            if re.fullmatch(r'^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', ncbi_id, re.IGNORECASE):
                add_cv_term_genes(ncbi_id, 'REFSEQ', gene, lab_strain)
                add_cv_term_genes(ncbi_id, 'NCBI', gene, lab_strain)
                name, locus = search_ncbi_for_gpr(ncbi_id)
                if (locus2id is not None) and (ncbi_id in locus2id.index):
                    locus = locus2id.loc[ncbi_id, 'LocusTag']

            # If identifier only contains numbers 
            # -> Get the corresponding data from the CarveMe input file
            elif re.fullmatch(r'^\d+$', ncbi_id, re.IGNORECASE):
                if id2locus_name is not None:
                    name, locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['name', 'locus_tag']].values[0]
                else: 
                    genes_missing_annotation.append(gene_id)
        
            # If identifier matches ncbiprotein ID pattern
            elif re.fullmatch(r'^(\w+\d+(\.\d+)?)|(NP_\d+)$', ncbi_id, re.IGNORECASE):
                add_cv_term_genes(ncbi_id, 'NCBI', gene, lab_strain)
                name, locus = search_ncbi_for_gpr(ncbi_id)
            
            # Catch all remaining cases that have no valid ID   
            else: 
                genes_missing_annotation.append(gene_id)
        
            # For lab strains use the locus tag from the annotation file   
            if lab_strain and id2locus_name is not None:
                locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['locus_tag']].values[0][0]
        
            if ncbi_id not in genes_missing_annotation:
                if name and locus:
                    gene.setName(name)
                    if len(locus) == 1: gene.setLabel(locus[0])
                    else: genes_missing_annotation.append(gene_id)
                else:
                    genes_missing_annotation.append(gene_id)
            
        gene.unsetNotes()
    if genes_missing_annotation:    
        genes_filename = f'{filename}_genes_missing_annotations_{str(date.today().strftime("%Y%m%d"))}.txt'
        logging.warning(f'''
                        For the following {len(genes_missing_annotation)} NCBI Protein IDs no annotation, name & label (locus tag) were found: {genes_missing_annotation}.
                        These IDs are saved to {genes_filename}
                        ''')
        with open(genes_filename, "w") as file:
             file.write(str(genes_missing_annotation))


#----------------------------  Functions to add additional URIs to GeneProducts ---------------------------------------#
def add_gp_id_from_gff(locus2id: pd.DataFrame, gene_list: list[GeneProduct]):
    """Adds URIs to GeneProducts based on locus tag to indentifier mapping

    Args:
        locus2id (pd.DataFrame): 
            Table mapping locus tags to their corresponding RefSeq identifiers
        gene_list (list[GeneProduct]): 
            libSBML ListOfGenes
    """
    locus2id.set_index('LocusTag')

    for gp in tqdm(gene_list):
        locus = gp.getLabel()

        if locus in locus2id.index:
            add_cv_term_genes(locus2id.loc[locus, 'ProteinID'][0].split(r'.')[0], 'REFSEQ', gp)
