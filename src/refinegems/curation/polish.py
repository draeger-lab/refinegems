#!/usr/bin/env python
"""General functions to polish a model

This module provides functionalities for an initial clean-up of models, including special functions for CarveMe models.

Since CarveMe version 1.5.1, the draft models from CarveMe contain pieces of information that are not correctly added to the 
annotations. To address this, this module includes the following functionalities:

    - Add URIs from the entity IDs to the annotation field for metabolites & reactions
    - Transfer URIs from the notes field to the annotations for metabolites & reactions
    - Add URIs from the GeneProduct IDs to the annotations

The functionalities for CarveMe models, along with the following further functionalities, are gathered in the main 
function :py:func:`~refinegems.curation.polish.polish`.

Further functionalities:

        - Setting boundary condition & constant for metabolites & reactions
        - Unit handling to add units & UnitDefinitions & to set units from parameters
        - Addition of default settings for compartments & metabolites
        - Addition of URIs to GeneProducts

            - via the RefSeq GFF from NCBI
            - via the KEGG API
            
        - Changing the CURIE pattern/CVTerm qualifier & qualifier type
        - Directionality control
"""

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel and Carolin Brune"
# @TODO Clean-up this module!
################################################################################
# requirements
################################################################################

import cobra
import logging
import pandas as pd
import re

from Bio import Entrez
from bioservices.kegg import KEGG
from colorama import init as colorama_init
from colorama import Fore
from datetime import date
from libsbml import Model as libModel
from libsbml import ListOfSpecies, ListOfReactions
from libsbml import GeneProduct, Species, Reaction, Unit, UnitDefinition
from libsbml import UNIT_KIND_MOLE, UNIT_KIND_GRAM, UNIT_KIND_LITRE, UNIT_KIND_SECOND
from libsbml import BQM_IS, BQM_IS_DERIVED_FROM, BQM_IS_DESCRIBED_BY

from tqdm.auto import tqdm
from typing import Union

from .miriam import polish_annotations, change_all_qualifiers
from ..utility.cvterms import add_cv_term_units, add_cv_term_metabolites, add_cv_term_reactions, add_cv_term_genes, metabol_db_dict, reaction_db_dict
from ..utility.db_access import search_ncbi_for_gpr
from ..utility.io import parse_gff_for_cds, parse_fasta_headers, load_a_table_from_database

################################################################################
# functions
################################################################################

# ------------------
# Polish - Main part
# ------------------
 
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
            for db in metabol_db_dict.keys():
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
            for db in reaction_db_dict.keys():
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


#--------------- Function to set boundary condition & constant for metabolites & reactions ----------------------------# 
def polish_metab_conditions(entity_list: Union[ListOfSpecies, ListOfReactions]):
    """Sets boundary condition and constant if not set for a metabolite

    Args:
        - entity_list (Union[ListOfSpecies, ListOfReactions]): 
            libSBML ListOfSpecies or ListOfReactions
    """
    match entity_list:
        case ListOfSpecies():
            for entity in entity_list:
                if not entity.getBoundaryCondition():
                    entity.setBoundaryCondition(False)
                if not entity.getConstant():
                    entity.setConstant(False)
        case ListOfReactions():
            pass
        case _:
            logging.warning(
                f'Unsupported entity_list type {type(entity_list)}. Please use only objects of type ListOfSpecies or ListOfReactions.'
                )


#------------------------------------ Functions to create units & UnitDefinitions -------------------------------------# 
def create_unit(
    model_specs: tuple[int], meta_id: str, kind: str, e: int, m: int, s: int, uri_is: str='', uri_idf: str=''
    ) -> Unit:
    """Creates unit for SBML model according to arguments

    Args:
        - model_specs (tuple):
            Level & Version of SBML model
        - meta_id (str): 
            Meta ID for unit (Neccessary for URI)
        - kind (str):
            Unit kind constant (see libSBML for available constants)
        - e (int): 
            Exponent of unit
        - m (int): 
            Multiplier of unit
        - s (int): 
            Scale of unit 
        - uri_is (str): 
            URI supporting the specified unit
        - uri_idf (str):
            URI supporting the derived from unit
      
    Returns:
        Unit: 
            libSBML unit object
    """
    unit = Unit(*model_specs)
    unit.setKind(kind)
    unit.setMetaId(f'meta_{meta_id}_{unit.getKind()}')
    unit.setExponent(e)
    unit.setMultiplier(m)
    unit.setScale(s)
    if uri_is:
        add_cv_term_units(uri_is, unit, BQM_IS)
    if uri_idf:
        add_cv_term_units(uri_idf, unit, BQM_IS_DERIVED_FROM)
    return unit


def create_unit_definition(model_specs: tuple[int], identifier: str, name: str, 
                           units: list[Unit]) -> UnitDefinition:
    """Creates unit definition for SBML model according to arguments
   
    Args:
        - model_specs (tuple): 
            Level & Version of SBML model
        - identifier (str):
            Identifier for the defined unit
        - name (str): 
            Full name of the defined unit
        - units (list): 
            All units the defined unit consists of
         
    Returns:
        UnitDefinition: 
            libSBML unit definition object
    """
    unit_definition = UnitDefinition(*model_specs)
    unit_definition.setId(identifier)
    unit_definition.setMetaId(f'meta_{identifier}')
    unit_definition.setName(name)
      
    # Iterate over all units provided for the unit definition & add the units
    for unit in units:
        unit_definition.addUnit(unit)
   
    return unit_definition


#----------------------------------- Function to create FBA unit & UnitDefinitions ------------------------------------#
def create_fba_units(model: libModel) -> list[UnitDefinition]:
    """Creates all fba units required for a constraint-based model
   
    Args:
        - model (libModel): 
            Model loaded with libSBML
         
    Returns:
        list: 
            List of libSBML UnitDefinitions
    """
    # Get model level & version for unit & unit definition
    model_specs = model.getLevel(), model.getVersion()
    
    # Create required units
    litre = create_unit(model_specs, 'litre', UNIT_KIND_LITRE, e=1, m=1, s=-3, uri_is='0000104', uri_idf='0000099')
    mole = create_unit(model_specs, 'mole_0', UNIT_KIND_MOLE, e=1, m=1, s=-3, uri_is='0000040', uri_idf='0000013')
    gram = create_unit(model_specs, 'per_gram_0', UNIT_KIND_GRAM, e=-1, m=1, s=0, uri_is='0000021')
    second = create_unit(model_specs, 'second_0', UNIT_KIND_SECOND, e=1, m=3600, s=0, uri_is='0000032', uri_idf='0000010')
    
    # Create unit definitions for hour & femto litre
    hour = create_unit_definition(
        model_specs, identifier='h', name='Hour', 
        units=[second]
    )
    femto_litre = create_unit_definition(
        model_specs, identifier='fL', name='Femto litres', 
        units=[litre]
    )
    
    # Create unit definitions for millimoles per gram dry weight (mmgdw) & mmgdw per hour
    mmgdw = create_unit_definition(
        model_specs, identifier='mmol_per_gDW', name='Millimoles per gram (dry weight)',
        units=[mole, gram]
    )

    # Create new units mole & gram to get new meta IDs 
    mole = create_unit(model_specs, 'mole_1', UNIT_KIND_MOLE, e=1, m=1, s=-3, uri_is='0000040', uri_idf='0000013')
    gram = create_unit(model_specs, 'per_gram_1', UNIT_KIND_GRAM, e=-1, m=1, s=0, uri_is='0000021')
    # Create new unit second to fit to per hour & get new meta ID
    second = create_unit(model_specs, 'second_1', UNIT_KIND_SECOND, e=-1, m=3600, s=0, uri_is='0000032', uri_idf='0000010')

    mmgdwh = create_unit_definition(
        model_specs, identifier='mmol_per_gDW_per_h', name='Millimoles per gram (dry weight) per hour',
        units=[mole, gram, second]
    )
    add_cv_term_units('pubmed:7986045', mmgdwh, BQM_IS_DESCRIBED_BY)

    return [mmgdwh, mmgdw, hour, femto_litre]


#--------------------------- Functions to print UnitDefinitions with the respective units -----------------------------#
def print_UnitDefinitions(contained_unit_defs: list[UnitDefinition]):
    """Prints a list of libSBML UnitDefinitions as XMLNodes
   
    Args:
        - contained_unit_defs (list): 
            List of libSBML UnitDefinition objects
    """
    for unit_def in contained_unit_defs:
        logging.info(unit_def.toXMLNode())


def print_remaining_UnitDefinitions(model: libModel, list_of_fba_units: list[UnitDefinition]):
    """Prints UnitDefinitions from the model that were removed as these were not contained in the list_of_fba_units

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - list_of_fba_units (list):  
            List of libSBML UnitDefinitions  
    """
       
    # Get all units already present in the model
    contained_unit_defs = [unit for unit in model.getListOfUnitDefinitions()]
         
    # Check if contained unit fits to one of the created fba units
    for unit_def in list_of_fba_units:
        for contained_unit_def in contained_unit_defs:
         
            current_id = contained_unit_def.getId()
               
            if UnitDefinition.areIdentical(unit_def, contained_unit_def):
                contained_unit_defs.remove(contained_unit_def)
                model.removeUnitDefinition(current_id)
   
    # Only print list if it contains UnitDefinitions         
    if contained_unit_defs:
        logging.info('''
        The following UnitDefinition objects were removed. 
        The reasoning is that
        \t(a) these UnitDefinitions are not contained in the UnitDefinition list of this program and
        \t(b) the UnitDefinitions defined within this program are handled as ground truth.
        Thus, the following UnitDefinitions are not seen as relevant for the model.
        ''')
        print_UnitDefinitions(contained_unit_defs)


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
            

#-------------------------- Functions to add default settings for compartments & metabolites --------------------------#            
def add_compartment_structure_specs(model: libModel):
    """| Adds the required specifications for the compartment structure
    | if not set (size & spatial dimension)
        
    Args:
        - model (libModel): 
            Model loaded with libSBML
    """ 
    for compartment in model.getListOfCompartments():
      
        if not compartment.isSetSize():
            compartment.setSize(float('NaN'))
         
        if not compartment.isSetSpatialDimensions():
            compartment.setSpatialDimensions(3)
         
        if any(
            (unit_id := re.fullmatch(r'fL', unit.getId(), re.IGNORECASE)) for unit in model.getListOfUnitDefinitions()
            ):
            if not (compartment.isSetUnits() and compartment.getUnits() == unit_id.group(0)):
                compartment.setUnits(unit_id.group(0))
         
         
def set_initial_amount(model: libModel):
    """Sets initial amount to all metabolites if not already set or if initial concentration is not set
    
    Args:
        - model (libModel): 
            Model loaded with libSBML
    """
    for species in model.getListOfSpecies():
      
      if not (species.isSetInitialAmount() or species.isSetInitialConcentration()):
         species.setInitialAmount(float('NaN'))


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

          
def add_gp_ids_from_KEGG(gene_list: list[GeneProduct], kegg_organism_id: str):
    """Adds KEGG gene & UniProt identifiers to the GeneProduct annotations

    Args:
        gene_list (list[GeneProduct]): 
            libSBML ListOfGenes
        kegg_organism_id (str): 
            Organism identifier in the KEGG database
    """
    k = KEGG()
    mapping_kegg_uniprot = k.conv('uniprot', kegg_organism_id)
    no_valid_kegg = []

    for gp in tqdm(gene_list):
    
        if gp.getId() != 'G_spontaneous':
            kegg_gene_id = f'{kegg_organism_id}:{gp.getLabel()}'
            
            try:
                uniprot_id = mapping_kegg_uniprot[kegg_gene_id]

                add_cv_term_genes(kegg_gene_id, 'KEGG', gp)
                add_cv_term_genes(uniprot_id.split(r'up:')[1], 'UNIPROT', gp)
                
            except KeyError:
                no_valid_kegg.append(gp.getLabel())
    
    if no_valid_kegg:      
        logging.info(f'The following {len(no_valid_kegg)} locus tags form no valid KEGG Gene ID: {no_valid_kegg}')

#--------------------------------------------------- Main function ----------------------------------------------------#
def polish(model: libModel, email: str, id_db: str, gff: str, protein_fasta: str, lab_strain: bool, 
           kegg_organism_id: str, path: str) -> libModel: 
    """Completes all steps to polish a model
    
    .. note:: So far only tested for models having either BiGG or VMH identifiers.

    Args:
        - model (libModel): 
            model loaded with libSBML
        - email (str): 
            E-mail for Entrez
        - id_db (str):
            Main database where identifiers in model come from
        - gff (str): 
            Path to RefSeq/Genbank GFF file of organism
        - protein_fasta (str): 
            File used as input for CarveMe
        - lab_strain (bool): 
            True if the strain was sequenced in a local lab
        - kegg_organism_id (str): 
            KEGG organism identifier
        - path (str): 
            Output path for incorrect annotations file(s)
    
    Returns:
        libModel: 
            Polished libSBML model
    """
    ### Set-up
    # Initialisation of colorama
    colorama_init(autoreset=True)
    
    # Filename for files tracking manual curation outcomes
    filename = f'{path}{model.getId()}'

    # Error/ invalid input handling
    if lab_strain and not protein_fasta:
        logging.error(Fore.LIGHTRED_EX + '''
                Setting the parameter lab_strain to True requires the provision of the protein FASTA file used as input for CarveMe.
                Otherwise, polish will not change anything for the GeneProducts.
                The header lines should look similar to the following line:
                >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
                It would also be a valid input if the header lines looked similar to the following line:
                >lcl|CP035291.1_prot_QCY37216.1_1 [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1]
                ''')
        return

    # Get ListOf objects
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    gene_list = model.getPlugin('fbc').getListOfGeneProducts()

    # Read GFF if provided
    if gff:
        locus2id = parse_gff_for_cds(gff, {'locus_tag':'LocusTag', 'protein_id':'ProteinID'})
        try: 
            locus2id = locus2id.explode('LocusTag').explode('ProteinID') # Replace (potentially single-entry) lists with their content
            locus2id.dropna(subset=['LocusTag'], axis=0, inplace=True) # If no locus tag exists, no mapping is possible
        except:
            mes = f'GFF does not contain the necessary information. Cannot be used.'
            logging.warning(mes)
            locus2id = None
    else: locus2id = None

    ### unit definition ###
    add_fba_units(model)
    set_default_units(model)
    set_units(model)
    add_compartment_structure_specs(model)
    set_initial_amount(model)
    
    ### improve metabolite, reaction and gene annotations ###
    add_metab(metab_list, id_db)
    add_reac(reac_list, id_db)
    cv_notes_metab(metab_list)
    cv_notes_reac(reac_list)
    # @DEBUG : comment out when fixing stuff in pollish_annotations (improves runtime) 
    cv_ncbiprotein(gene_list, email, locus2id, protein_fasta, filename, lab_strain) 
    

    ### add additional URIs to GeneProducts ###
    if locus2id is not None: add_gp_id_from_gff(locus2id, gene_list)
    if kegg_organism_id: add_gp_ids_from_KEGG(gene_list, kegg_organism_id)
    
    ### set boundaries and constants ###
    polish_metab_conditions(metab_list)
    polish_metab_conditions(reac_list)
    
    ### MIRIAM compliance of CVTerms ###
    print('Remove duplicates & transform all CURIEs to the new identifiers.org pattern (: between db and ID):')
    model = polish_annotations(model, True, filename)
    print('Changing all qualifiers to be MIRIAM compliant:')
    model = change_all_qualifiers(model, lab_strain)
    
    return model



# ----------------------
# Directionality Control
# ----------------------

def check_direction(model:cobra.Model,data:Union[pd.DataFrame,str]) -> cobra.Model:
    """Check the direction of reactions by searching for matching MetaCyc,
    KEGG and MetaNetX IDs as well as EC number in a downloaded BioCyc (MetaCyc)
    database table or dataFrame (need to contain at least the following columns:
    Reactions (MetaCyc ID),EC-Number,KEGG reaction,METANETX,Reaction-Direction.

    Args:
        model (cobra.Model): 
            The model loaded with COBRApy.
        data (pd.DataFrame | str): 
            Either a pandas DataFrame or a path to a CSV file
            containing the BioCyc smart table.

    Raises:
        - TypeError: Unknown data type for parameter data

    Returns:
        cobra.Model: 
            The edited model.
    """
    
    match data:
        # already a DataFrame
        case pd.DataFrame():
            pass
        case str():
            # load from a table
            data = pd.read_csv(data, sep='\t')
            # rewrite the columns into a better comparable/searchable format
            data['KEGG reaction'] = data['KEGG reaction'].str.extract(r'.*>(R\d*)<.*')
            data['METANETX']      = data['METANETX'].str.extract(r'.*>(MNXR\d*)<.*')
            data['EC-Number']     = data['EC-Number'].str.extract(r'EC-(.*)')
        case _:
            mes = f'Unknown data type for parameter data: {type(data)}'
            raise TypeError(mes)

    # check direction
    # --------------------
    for r in model.reactions:

        direction = None
        # easy case: metacyc is already (corretly) annotated
        if 'metacyc.reaction' in r.annotation and len(data[data['Reactions'] == r.annotation['metacyc.reaction']]) != 0:
            direction = data[data['Reactions'] == r.annotation['metacyc.reaction']]['Reaction-Direction'].iloc[0]
            r.notes['BioCyc direction check'] = F'found {direction}'
        # complicated case: no metacyc annotation
        else:
            annotations = []

            # collect matches
            if 'kegg.reaction' in r.annotation and r.annotation['kegg.reaction'] in data['KEGG reaction'].tolist():
                annotations.append(data[data['KEGG reaction'] == r.annotation['kegg.reaction']]['Reactions'].tolist())
            if 'metanetx.reaction' in r.annotation and r.annotation['metanetx.reaction'] in data['METANETX'].tolist():
                annotations.append(data[data['METANETX'] == r.annotation['metanetx.reaction']]['Reactions'].tolist())
            if 'ec-code' in r.annotation and r.annotation['ec-code'] in data['EC-Number'].tolist():
                annotations.append(data[data['EC-Number'] == r.annotation['ec-code']]['Reactions'].tolist())

            # check results
            # no matches
            if len(annotations) == 0:
                r.notes['BioCyc direction check'] = 'not found'

            # matches found
            else:
                # built intersection
                intersec = set(annotations[0]).intersection(*annotations)
                # case 1: exactly one match remains
                if len(intersec) == 1:
                    entry = intersec.pop()
                    direction = data[data['Reactions'] == entry]['Reaction-Direction'].iloc[0]
                    r.annotation['metacyc.reaction'] = entry
                    r.notes['BioCyc direction check'] = F'found {direction}'

                # case 2: multiple matches found -> inconclusive
                else:
                    r.notes['BioCyc direction check'] = F'found, but inconclusive'

        # update direction if possible and needed
        if not pd.isnull(direction):
            if 'REVERSIBLE' in direction:
                # set reaction as reversible by setting default values for upper and lower bounds
                r.lower_bound = cobra.Configuration().lower_bound
            elif 'RIGHT-TO-LEFT' in direction:
                # invert the default values for the boundaries
                r.lower_bound = cobra.Configuration().lower_bound
                r.upper_bound = 0.0
            elif 'LEFT-To-RIGHT' in direction:
                # In case direction was already wrong
                r.lower_bound = 0.0
                r.upper_bound = cobra.Configuration().upper_bound
            else:
                # left to right case is the standard for adding reactions
                # = nothing left to do
                continue

    return model
