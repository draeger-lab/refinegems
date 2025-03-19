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
from libsbml import Model as libModel
from libsbml import GeneProduct, Species, Reaction

from tqdm.auto import tqdm

from ..utility.cvterms import add_cv_term_metabolites, add_cv_term_reactions, add_cv_term_genes, DB2PREFIX_METABS, DB2PREFIX_REACS
from ..utility.db_access import search_ncbi_for_gpr
from ..utility.entities import create_fba_units, print_remaining_UnitDefinitions
from ..utility.io import parse_fasta_headers, load_a_table_from_database

################################################################################
# functions
################################################################################
 
#----------- Functions to add URIs from the entity IDs to the annotation field for metabolites & reactions ------------#       
def add_metab(entity_list: list[Species], id_db: str):
    """| Adds the ID of metabolites as URI to the annotation field
    | For a VMH model, additionally, the corresponding BiGG IDs are added! 
	
    (Currently, only BiGG & VMH IDs supported!)

    Args:
        - entity_list (list): 
            libSBML ListOfSpecies
        - id_db (str): 
            Name of the database of the IDs contained in a model.
    """
    # @TODO Use cobra.io.sbml._f_specie
    vmh_cut_pattern = r'__\d+_' # To extract BiGG identifier
    
    if id_db == 'VMH':
        bigg_metabs_ids = load_a_table_from_database('SELECT universal_bigg_id FROM bigg_metabolites')
        bigg_metabs_ids = bigg_metabs_ids['universal_bigg_id'].tolist()

    for entity in entity_list:
        
        # Get ID and remove 'M_'
        current_id = entity.getId()
        # @TODO Use here cobra.io.sbml._f_specie, also removes prefix 'M_'
        current_id = current_id[2:]
        
        # Unset annotations if no CV terms exist
        if entity.getNumCVTerms() == 0:
            entity.unsetAnnotation()

        # Use current_id as metaid if no metaid is present   
        if not entity.isSetMetaId():
            entity.setMetaId(f'meta_M_{current_id}')
        
        # If database 'VMH' was specified, extract BiGG ID with vmh_cut_pattern
        if id_db == 'VMH':
            # id_for_anno = current_id[:-3]
            current_id_cut = re.split(vmh_cut_pattern, current_id)
            current_id = ''.join(current_id_cut[:2])
        
        # Remove compartment specifier from ID for annotation   
        id_for_anno = current_id[:-2]

		# Add ID as URI to annotation
        add_cv_term_metabolites(id_for_anno, id_db, entity)

        # Add BiGG ID to annotation, additionally, if valid BiGG ID
        if (id_db == 'VMH') and (id_for_anno in bigg_metabs_ids):
                add_cv_term_metabolites(id_for_anno, 'BIGG', entity)
            
           
def add_reac(entity_list: list[Reaction], id_db: str):
    """| Adds the ID of reactions as URI to the annotation field
       | (Currently, only BiGG & VMH IDs supported!)

        Args:
            - entity_list (list): 
                libSBML ListOfReactions
            - id_db (str): 
                Name of the database of the IDs contained in a model             
    """
    # @TODO Use cobra.io.sbml._f_reaction
    if id_db == 'VMH':
       bigg_reacs_ids = load_a_table_from_database('SELECT bigg_id FROM bigg_reactions')
       bigg_reacs_ids = bigg_reacs_ids['bigg_id'].tolist()

    # Use regex to generalise check for growth/biomass reaction
    regex = r'growth|_*biomass\d*_*'
    
    for entity in entity_list:
        
        # Get ID and remove 'R_'
        current_id = entity.getId()
        if current_id[:2] == 'R_':
            current_id = current_id[2:]
        
        # Readjusted to fit to the metaid pattern from the metabolites      
        if not entity.isSetMetaId():
            entity.setMetaId(f'meta_R_{current_id}')

        if not re.fullmatch(regex, current_id, re.IGNORECASE):
            
            # Unset annotations if no CV terms exist
            if entity.getNumCVTerms() == 0:
                entity.unsetAnnotation()
            
            # Add ID as URI to annotation   
            add_cv_term_reactions(current_id, id_db, entity)

            # Add BiGG ID to annotation, additionally, if valid BiGG ID
            if (id_db == 'VMH') and (current_id in bigg_reacs_ids):
                add_cv_term_reactions(current_id, 'BIGG', entity)


#----------- Functions to transfer URIs from the notes field to the annotations for metabolites & reactions -----------# 
def cv_notes_metab(species_list: list[Species]):
    """Checks the notes field for information which should be in the annotation field.
    Removes entry from notes and adds it as URL to the CVTerms of a metabolite

    Args:
        - species_list (list): 
            libSBML ListOfSpecies
    """

    for species in species_list:
        if not species.isSetMetaId():
            species.setMetaId('meta_' + species.getId())
        notes_list = []
        elem_used = []
        notes_string = species.getNotesString().split(r'\n')
        for elem in notes_string:
            for db in DB2PREFIX_METABS.keys():
                if '<p>' + db in elem:
                    elem_used.append(elem)
                    # @DEBUG print(elem.strip()[:-4].split(r': ')[1])
                    fill_in = re.split(r':\s*', elem.strip()[:-4])[1]
                    if (';') in fill_in and not re.search(r'inchi', db, re.IGNORECASE):
                        entries = fill_in.split(r';')
                        for entry in entries:
                            if not re.fullmatch(r'^nan$', entry.strip(), re.IGNORECASE):
                                add_cv_term_metabolites(entry.strip(), db, species)
                    else:
                        if not re.fullmatch(r'^nan$', fill_in, re.IGNORECASE):
                            add_cv_term_metabolites(fill_in, db, species)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        species.unsetNotes()
        species.setNotes(new_notes)
        # @DEBUG print(species.getAnnotationString())


def cv_notes_reac(reaction_list: list[Reaction]):
    """Checks the notes field for information which should be in the annotation field.
    Removes entry from notes and adds it as URL to the CVTerms of a reaction

    Args:
        - reaction_list (list): 
            libSBML ListOfReactions
    """

    for reaction in reaction_list:
        if not reaction.isSetMetaId():
            reaction.setMetaId('meta_' + reaction.getId())
        notes_list = []
        elem_used = []
        notes_string = reaction.getNotesString().split(r'\n')

        for elem in notes_string:
            for db in DB2PREFIX_REACS.keys():
                if '<p>' + db in elem:
                    elem_used.append(elem)
                    # @DEBUG print(elem.strip()[:-4].split(r': ')[1])
                    fill_in = re.split(r':\s*', elem.strip()[:-4])[1]
                    if (';') in fill_in:
                        entries = fill_in.split(r';')
                        for entry in entries:
                            if not re.fullmatch(r'^nan$', entry.strip(), re.IGNORECASE):
                                add_cv_term_reactions(entry.strip(), db, reaction)
                    else:
                        if not re.fullmatch(r'^nan$', fill_in, re.IGNORECASE):
                            add_cv_term_reactions(fill_in, db, reaction)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        reaction.unsetNotes()
        reaction.setNotes(new_notes)


#-------------------------------------- Function to add units & UnitDefinitions ---------------------------------------#
# @TODO Redo unit function setup
# @IDEA
# validate_fba_units
# correct_fba_units -> Only if necessary
# set_fba_units
def add_fba_units(model: libModel):
    """Adds:
    
       - mmol per gDW per h
       - mmol per gDW 
       - hour (h)
       - femto litre (fL)
       
       to the list of unit definitions (needed for FBA)

    Args:
        - model (libModel): 
            Model loaded with libSBML
    """
    list_of_fba_units = create_fba_units(model)
    
    # If list of unit definitions is not empty, replace all units identical to list_of_fba_units
    # & Print the remaining unit definitions
    if model.getListOfUnitDefinitions():
        print_remaining_UnitDefinitions(model, list_of_fba_units)
    
    for unit_def in list_of_fba_units:
        model.getListOfUnitDefinitions().append(unit_def)


def set_default_units(model: libModel):
    """Sets default units of model

    Args:
        - model (libModel): 
            Model loaded with libSBML
    """ 
    for unit in model.getListOfUnitDefinitions():
      
        unit_id = unit.getId()
      
        if re.fullmatch(r'mmol_per_gDW', unit_id, re.IGNORECASE):
         
            if not (model.isSetExtentUnits() and model.getExtentUnits() == unit_id):
                model.setExtentUnits(unit_id)
            
            if not (model.isSetSubstanceUnits() and model.getSubstanceUnits() == unit_id):
                model.setSubstanceUnits(unit_id)
         
        if not (model.isSetTimeUnits() and model.getTimeUnits() == unit_id) and re.fullmatch(r'hr?', unit_id, re.IGNORECASE):
            model.setTimeUnits(unit_id)
         
        if not (model.isSetVolumeUnits() and model.getVolumeUnits() == unit_id) and re.fullmatch(r'fL', unit_id, re.IGNORECASE):
            model.setVolumeUnits(unit_id)


#--------------------------------------- Function to set units from parameters ----------------------------------------#
def set_units(model: libModel):
    """Sets units of parameters in model

    Args:
        - model (libModel): 
            Model loaded with libSBML
    """
    for param in model.getListOfParameters(): # needs to be added to list of unit definitions aswell
        if any(
            (unit_id := re.fullmatch(r'mmol_per_gDW_per_hr?', unit.getId(), re.IGNORECASE)) 
            for unit in model.getListOfUnitDefinitions()
            ):
            if not (param.isSetUnits() and param.getUnits() == unit_id.group(0)):
                param.setUnits(unit_id.group(0))


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
