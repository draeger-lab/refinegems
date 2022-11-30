#!/usr/bin/env python
""" Can be used to polish a model (created with CarveMe v.1.5.1)

The newer version of CarveMe leads to some irritations in the model, these scripts enable for example the addition of BiGG Ids to the annotations as well as a correct formatting of the annotations.
"""

import re
from libsbml import *
from Bio import Entrez, SeqIO
from tqdm.auto import tqdm
from refinegems.cvterms import add_cv_term_units, add_cv_term_metabolites, add_cv_term_reactions, add_cv_term_genes, metabol_db_dict, reaction_db_dict
from refinegems.load import write_to_file

__author__ = "Famke Baeuerle"
__author__ = 'Gwendolyn O. Gusak'
        
        
def add_metab(entity_list, id_db: str):
    """adds the ID of metabolites as URI to the annotation field
       For a VMH model, additionally, the corresponding BiGG IDs are added! 
		(Currently, only BiGG & VMH IDs supported!)

    Args:
        entity_list (list): libSBML ListOfSpecies
    """
    vmh_cut_pattern = '__\d+_' # To extract BiGG identifier
    
    for entity in entity_list:
        
        # Get ID and remove 'M_'
        current_id = entity.getId()
        current_id = current_id[2:]
        
        # Unset annotations if no CV terms exist
        if entity.getNumCVTerms() == 0:
           entity.unsetAnnotation()
        
        # If database 'VMH' was specified, extract BiGG ID with vmh_cut_pattern
        if id_db == 'VMH':
           current_id_cut = re.split(vmh_cut_pattern, current_id)
           current_id = ''.join(current_id_cut[:2])
           
        # Use current_id as metaid if no metaid is present   
        if not entity.isSetMetaId():
            entity.setMetaId(f'meta_M_{current_id}')
        
        # Remove compartment specifier from ID for annotation   
        id_for_anno = current_id[:-2]

		  # Add ID as URI to annotation
        add_cv_term_metabolites(id_for_anno, id_db, entity)
        
        if id_db == 'VMH':
            # Add BiGG ID to annotation, additionally
            add_cv_term_metabolites(id_for_anno, 'BIGG', entity)
            
            
def add_reac(entity_list, id_db: str):
   """adds the ID of reactions as URI to the annotation field
		(Currently, only BiGG & VMH IDs supported!)

	Args:
		entity_list (list): libSBML ListOfReactions
	"""
   # Use regex to generalise check for growth/biomass reaction
   regex = 'growth|_*biomass\d*_*'
   
   for entity in entity_list:
      
      # Get ID and remove 'R_'
      current_id = entity.getId()
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


def cv_notes_metab(species_list):
    """checks the notes field for information which should be in the annotation field
       removes entry from notes and adds it as URL to the CVTerms of a metabolite

    Args:
        species_list (list): libSBML ListOfSpecies
    """

    for species in species_list:
        if not species.isSetMetaId():
            species.setMetaId('meta_' + species.getId())
        notes_list = []
        elem_used = []
        notes_string = species.getNotesString().split('\n')
        for elem in notes_string:
            for db in metabol_db_dict.keys():
                if '<p>' + db in elem:
                    elem_used.append(elem)
                    #print(elem.strip()[:-4].split(': ')[1])
                    fill_in = re.split(':\s*', elem.strip()[:-4])[1]
                    if (';') in fill_in and db != 'INCHI':
                        entries = fill_in.split(';')
                        for entry in entries:
                            add_cv_term_metabolites(entry.strip(), db, species)
                    else:
                        add_cv_term_metabolites(fill_in, db, species)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        species.unsetNotes()
        species.setNotes(new_notes)
        #print(species.getAnnotationString())


def cv_notes_reac(reaction_list):
    """checks the notes field for information which should be in the annotation field
       removes entry from notes and adds it as URL to the CVTerms of a reaction

    Args:
        reaction_list (list): libSBML ListOfReactions
    """

    for reaction in reaction_list:
        if not reaction.isSetMetaId():
            reaction.setMetaId('meta_' + reaction.getId())
        notes_list = []
        elem_used = []
        notes_string = reaction.getNotesString().split('\n')

        for elem in notes_string:
            for db in reaction_db_dict.keys():
                if '<p>' + db in elem:
                    elem_used.append(elem)
                    #print(elem.strip()[:-4].split(': ')[1])
                    fill_in = re.split(':\s*', elem.strip()[:-4])[1]
                    if (';') in fill_in:
                        entries = fill_in.split(';')
                        for entry in entries:
                            add_cv_term_reactions(entry.strip(), db, reaction)
                    else:
                        add_cv_term_reactions(fill_in, db, reaction)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        reaction.unsetNotes()
        reaction.setNotes(new_notes)


def polish_entities(entity_list, metabolite):
    """sets boundary condition and constant if not set for a metabolite

    Args:
        entity_list (list): libSBML ListOfSpecies or ListOfReactions
        metabolite (boolean): flag to determine whether entity = metabolite
    """
    for entity in entity_list:
        if metabolite:  # some polishing
            if not entity.getBoundaryCondition():
                entity.setBoundaryCondition(False)
            if not entity.getConstant():
                entity.setConstant(False) 
 
                
def create_unit(model_specs: tuple[int], meta_id: str, kind: str, e: int, m: int, s: int, uri_is: str='', uri_idf: str='') -> Unit:
   '''Creates unit for SBML model according to parameters

      Params:
         model_specs (tuple):    Level & Version of SBML model
         meta_id (str):          Meta ID for unit (Neccessary for URI)
         kind (str):             Unit kind constant (see libSBML for available constants)
         e (int):                Exponent of unit
         m (int):                Multiplier of unit
         s (int):                Scale of unit 
         uri_is (str):           URI supporting the specified unit
         uri_idf (str):          URI supporting the derived from unit
      
      Return:
         libSBML unit object
   '''
   unit = Unit(*model_specs)
   unit.setKind(kind)
   unit.setMetaId(f'meta_{meta_id}_{unit.getKind()}')
   unit.setExponent(e)
   unit.setMultiplier(m)
   unit.setScale(s)
   if uri_is:
      add_cv_term_units(uri_is, unit)
   if uri_idf:
      add_cv_term_units(uri_idf, unit, True)
   return unit


def create_unit_definition(model_specs: tuple[int], identifier: str, name: str, 
                           units: list[Unit]) -> UnitDefinition:
   '''Creates unit definition for SBML model according to parameters
   
      Params:
         model_specs (tuple): Level & Version of SBML model
         identifier (str):    Identifier for the defined unit
         name (str):          Full name of the defined unit
         units (list):        All units the defined unit consists of
         
      Return:
         libSBML unit definition object
   '''
   unit_definition = UnitDefinition(*model_specs)
   unit_definition.setId(identifier)
   unit_definition.setMetaId(f'meta_{identifier}')
   unit_definition.setName(name)
      
   # Iterate over all units provided for the unit definition & add the units
   for unit in units:
      unit_definition.addUnit(unit)
   
   return unit_definition


def create_fba_units(model: Model) -> list[UnitDefinition]:
   '''Creates all fba units required for a constraint-based model
   
      Params:
         model (Model): SBML model loaded with libSBML
         
      Return:
         list of libSBML UnitDefinitions
   '''
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
   
   return [hour, femto_litre, mmgdw, mmgdwh]


def get_missing_fba_units(model: Model, list_of_fba_units: list[UnitDefinition]) -> list[UnitDefinition]:
   '''Retrieves a list of the missing units in a model 

      Params:
         model (Model): SBML model loaded with libSBML
         list_of_fba_units (list):  list of libSBML UnitDefinitions
         
      Return:
         list of libSBML UnitDefinition    
   '''
   polished_mmgdwh = False
   model_specs = model.getLevel(), model.getVersion()
   id_mapper = {'h': 'hr?', 'fL': 'fL', 'mmol_per_gDW': 'mmol_per_gDW', 'mmol_per_gDW_per_h': 'mmol_per_gDW_per_hr?'}
   new_list = []
       
   # Get all units already present in the model
   contained_unit_defs = [unit.getId() for unit in model.getListOfUnitDefinitions()]
         
   # Check if contained unit fits to one of the created fba units
   for unit_def in list_of_fba_units:
      
      is_matched = False
      
      for contained_unit_def in contained_unit_defs:
         
         # Remove units added by ModelPolisher & Replace with the ones from the list
         if 'substance' == contained_unit_def:
            contained_unit_defs.remove(contained_unit_def)
            model.removeUnitDefinition(contained_unit_def)
         
         if 'time' == contained_unit_def:
            contained_unit_defs.remove(contained_unit_def)
            model.removeUnitDefinition(contained_unit_def)
            
         if re.fullmatch('mmol_per_gDW_per_hr?', contained_unit_def, re.IGNORECASE) and not polished_mmgdwh:
            polished_mmgdwh = True
            mmgdwh = model.getUnitDefinition(contained_unit_def)
            
            if not mmgdwh.isSetMetaId():
               mmgdwh.setMetaId(f'meta_{mmgdwh.getId()}')
            
            if not mmgdwh.isSetName():
               mmgdwh.setName('Millimoles per gram (dry weight) per hour')
            
            # Adjust unit list
            mmgdwh.getListOfUnits().clear()
            mmgdwh.addUnit(create_unit(model_specs, 'mole_1', UNIT_KIND_MOLE, e=1, m=1, s=-3, uri_is='0000040', uri_idf='0000013'))
            mmgdwh.addUnit(create_unit(model_specs, 'per_gram_1', UNIT_KIND_GRAM, e=-1, m=1, s=0, uri_is='0000021'))
            mmgdwh.addUnit(create_unit(model_specs, 'second_1', UNIT_KIND_SECOND, e=-1, m=3600, s=0, uri_is='0000032', uri_idf='0000010'))
               
         if re.fullmatch(id_mapper.get(unit_def.getId()), contained_unit_def, re.IGNORECASE):
            is_matched = True       
            
      if not is_matched and (unit_def not in new_list):
         new_list.append(unit_def)
               
   return new_list


def add_fba_units(model: Model):
    """adds 
         - hour (h)
         - femto litre (fL)
         - mmol per gDW
         - mmol per gDW per h 
       to the list of unit definitions (needed for FBA)

    Args:
        model (libsbml-model): model loaded with libsbml
    """
    list_of_fba_units = create_fba_units(model)
    
    # If list of unit definitions is not empty, only add missing units
    if model.getListOfUnitDefinitions():
       list_of_fba_units = get_missing_fba_units(model, list_of_fba_units)
    
    for unit_def in list_of_fba_units:
         model.getListOfUnitDefinitions().append(unit_def)
          

def set_default_units(model: Model):
   '''Sets default units of model

      Params:
         model (Model): SBML model loaded with libSBML
   '''
   for unit in model.getListOfUnitDefinitions():
      
      unit_id = unit.getId()
      
      if re.fullmatch('mmol_per_gDW', unit_id, re.IGNORECASE):
         
         if not (model.isSetExtentUnits() and model.getExtentUnits() == unit_id):
            model.setExtentUnits(unit_id)
            
         if not (model.isSetSubstanceUnits() and model.getSubstanceUnits() == unit_id):
            model.setSubstanceUnits(unit_id)
         
      if not (model.isSetTimeUnits() and model.getTimeUnits() == unit_id) and re.fullmatch('hr?', unit_id, re.IGNORECASE):
         model.setTimeUnits(unit_id)
         
      if not (model.isSetVolumeUnits() and model.getVolumeUnits() == unit_id) and re.fullmatch('fL', unit_id, re.IGNORECASE):
         model.setVolumeUnits(unit_id)


def set_units(model):
    """Sets units of parameters in model

    Args:
        model (libsbml-model): model loaded with libsbml
    """
    for param in model.getListOfParameters(): # needs to be added to list of unit definitions aswell
        if param.isSetUnits() == False:
           if any((unit_id := re.fullmatch('mmol_per_gDW_per_hr?', unit.getId(), re.IGNORECASE)) for unit in model.getListOfUnitDefinitions()):
              param.setUnits(unit_id.group(0))
            
            
def add_compartment_structure_specs(model: Model):
   ''' Adds the required specifications for the compartment structure
       if not set (size & spatial dimension)
   '''
   for compartment in model.getListOfCompartments():
      
      if not compartment.isSetSize():
         compartment.setSize(float('NaN'))
         
      if not compartment.isSetSpatialDimensions():
         compartment.setSpatialDimensions(3)
         
         
def set_initial_amount(model: Model):
   
   for species in model.getListOfSpecies():
      
      if not (species.isSetInitialAmount() or species.isSetInitialConcentration()):
         species.setInitialAmount(float('NaN'))


def cv_ncbiprotein(gene_list, email):
    """Adds NCBI Id to genes as annotation

    Args:
        gene_list (list): libSBML ListOfGenes
        email (string): User Email to access the Entrez database
    """
    Entrez.email = email

    def get_name_locus_tag(ncbi_id):
        handle = Entrez.efetch(
            db="protein",
            id=ncbi_id,
            rettype="gbwithparts",
            retmode='text')
        records = SeqIO.parse(handle, "gb")

        for i, record in enumerate(records):
            if (ncbi_id[0] == 'W'):
                return record.description, ncbi_id
            else:
                for feature in record.features:
                    if feature.type == "CDS":
                        return record.description, feature.qualifiers["locus_tag"][0]

    print('Setting CVTerms and removing notes for all genes:')
    for gene in tqdm(gene_list):
        if not gene.isSetMetaId():
            gene.setMetaId('meta_' + gene.getId())
        
        if (gene.getId()[2] == 'W'): #addition to work with KC-Na-01
            entry = gene.getId()[2:-2]
            add_cv_term_genes(entry, 'NCBI', gene)
            name, locus = get_name_locus_tag(entry)
            gene.setName(name)
            gene.setLabel(locus)
        
        else:
            id_string = gene.getId().split('_')
            for entry in id_string:
                if (entry[:1] == 'Q'):  # maybe change this to user input
                    add_cv_term_genes(entry, 'NCBI', gene)
                    name, locus = get_name_locus_tag(entry)
                    gene.setName(name)
                    gene.setLabel(locus)

        gene.unsetNotes()
        
def polish(model, new_filename, email, id_db):
    """completes all steps to polish a model
         (Tested for models having either BiGG or VMH identifiers.)

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    gene_list = model.getPlugin('fbc').getListOfGeneProducts()

    add_fba_units(model)
    set_default_units(model)
    set_units(model)
    add_compartment_structure_specs(model)
    set_initial_amount(model)
    add_metab(metab_list, id_db)
    add_reac(reac_list, id_db)
    cv_notes_metab(metab_list)
    cv_notes_reac(reac_list)
    cv_ncbiprotein(gene_list, email)
    polish_entities(metab_list, metabolite=True)
    polish_entities(reac_list, metabolite=False)

    write_to_file(model, new_filename)
