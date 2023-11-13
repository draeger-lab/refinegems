#!/usr/bin/env python
"""Can be used to polish a model (created with CarveMe v.1.5.1)

The newer version of CarveMe leads to some irritations in the model, these scripts enable for example the addition of BiGG Ids to the annotations as well as a correct formatting of the annotations.
"""

import re, logging
from bioregistry import is_valid_identifier, manager, normalize_parsed_curie, get_identifiers_org_iri
from libsbml import Model as libModel
from libsbml import Species, Reaction, Unit, UnitDefinition, SBase, UNIT_KIND_MOLE, UNIT_KIND_GRAM, UNIT_KIND_LITRE, UNIT_KIND_SECOND, MODEL_QUALIFIER, BQM_IS, BQM_IS_DERIVED_FROM, BQM_IS_DESCRIBED_BY, BIOLOGICAL_QUALIFIER, BQB_IS, BQB_HAS_PROPERTY, BQB_IS_HOMOLOG_TO, BiolQualifierType_toString, ModelQualifierType_toString
from Bio import Entrez
from tqdm.auto import tqdm
from functools import reduce
from sortedcontainers import SortedDict, SortedSet
from refinegems.cvterms import add_cv_term_units, add_cv_term_metabolites, add_cv_term_reactions, add_cv_term_genes, generate_cvterm, metabol_db_dict, reaction_db_dict, MIRIAM, OLD_MIRIAM
from refinegems.io import search_ncbi_for_gpr, parse_fasta_headers, parse_dict_to_dataframe, load_a_table_from_database
from colorama import init as colorama_init
from colorama import Fore, Style
from datetime import date

__author__ = "Famke Baeuerle and Gwendolyn O. Döbel"


#------------------------------------------------ Constant variables --------------------------------------------------#
BIOCYC_TIER1_DATABASES_PREFIXES = ['META', 'ECO', 'ECOLI', 'HUMAN']
    
 
#----------- Functions to add URIs from the entity IDs to the annotation field for metabolites & reactions ------------#       
def add_metab(entity_list: list[Species], id_db: str):
    """| Adds the ID of metabolites as URI to the annotation field
       | For a VMH model, additionally, the corresponding BiGG IDs are added! 
	   | (Currently, only BiGG & VMH IDs supported!)

    Args:
        - entity_list (list): libSBML ListOfSpecies
        - id_db (str): Name of the database of the IDs contained in a model 
    """
    vmh_cut_pattern = '__\d+_' # To extract BiGG identifier
    
    if id_db == 'VMH':
        bigg_metabs_ids = load_a_table_from_database('SELECT universal_bigg_id FROM bigg_metabolites')
        bigg_metabs_ids = bigg_metabs_ids['universal_bigg_id'].tolist()

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
            # Check if valid BiGG ID
            if id_for_anno in bigg_metabs_ids:
                # Add BiGG ID to annotation, additionally
                add_cv_term_metabolites(id_for_anno, 'BIGG', entity)
            
           
def add_reac(entity_list: list[Reaction], id_db: str):
    """| Adds the ID of reactions as URI to the annotation field
       | (Currently, only BiGG & VMH IDs supported!)

        Args:
            - entity_list (list): libSBML ListOfReactions
            - id_db (str): Name of the database of the IDs contained in a model             
    """
    if id_db == 'VMH':
       bigg_reacs_ids = load_a_table_from_database('SELECT bigg_id FROM bigg_reactions')
       bigg_reacs_ids = bigg_reacs_ids['bigg_id'].tolist()

    # Use regex to generalise check for growth/biomass reaction
    regex = 'growth|_*biomass\d*_*'
    
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

            # If VMH ID == BiGG ID, add BiGG ID as well
            if (id_db == 'VMH') and (current_id in bigg_reacs_ids):
                add_cv_term_reactions(current_id, 'BIGG', entity)


#----------- Functions to transfer URIs from the notes field to the annotations for metabolites & reactions -----------# 
def cv_notes_metab(species_list: list[Species]):
    """| Checks the notes field for information which should be in the annotation field
       | removes entry from notes and adds it as URL to the CVTerms of a metabolite

    Args:
        - species_list (list): libSBML ListOfSpecies
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
                    if (';') in fill_in and not re.search('inchi', db, re.IGNORECASE):
                        entries = fill_in.split(';')
                        for entry in entries:
                            if not re.fullmatch('^nan$', entry.strip(), re.IGNORECASE):
                                add_cv_term_metabolites(entry.strip(), db, species)
                    else:
                        if not re.fullmatch('^nan$', fill_in, re.IGNORECASE):
                            add_cv_term_metabolites(fill_in, db, species)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        species.unsetNotes()
        species.setNotes(new_notes)
        #print(species.getAnnotationString())


def cv_notes_reac(reaction_list: list[Reaction]):
    """| Checks the notes field for information which should be in the annotation field
       | removes entry from notes and adds it as URL to the CVTerms of a reaction

    Args:
        - reaction_list (list): libSBML ListOfReactions
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
                            if not re.fullmatch('^nan$', entry.strip(), re.IGNORECASE):
                                add_cv_term_reactions(entry.strip(), db, reaction)
                    else:
                        if not re.fullmatch('^nan$', fill_in, re.IGNORECASE):
                            add_cv_term_reactions(fill_in, db, reaction)

        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)

        # Adding new, shortened notes
        new_notes = ' '.join([str(elem) + '\n' for elem in notes_list])
        reaction.unsetNotes()
        reaction.setNotes(new_notes)


#--------------- Function to set boundary condition & constant for metabolites & reactions ----------------------------# 
def polish_entities(entity_list: list, metabolite: bool):
    """Sets boundary condition and constant if not set for a metabolite

    Args:
        - entity_list (list): libSBML ListOfSpecies or ListOfReactions
        - metabolite (boolean): flag to determine whether entity = metabolite
    """
    for entity in entity_list:
        if metabolite:  # some polishing
            if not entity.getBoundaryCondition():
                entity.setBoundaryCondition(False)
            if not entity.getConstant():
                entity.setConstant(False) 


#------------------------------------ Functions to create units & UnitDefinitions -------------------------------------# 
def create_unit(
    model_specs: tuple[int], meta_id: str, kind: str, e: int, m: int, s: int, uri_is: str='', uri_idf: str=''
    ) -> Unit:
    """Creates unit for SBML model according to arguments

    Args:
        - model_specs (tuple):    Level & Version of SBML model
        - meta_id (str):          Meta ID for unit (Neccessary for URI)
        - kind (str):             Unit kind constant (see libSBML for available constants)
        - e (int):                Exponent of unit
        - m (int):                Multiplier of unit
        - s (int):                Scale of unit 
        - uri_is (str):           URI supporting the specified unit
        - uri_idf (str):          URI supporting the derived from unit
      
    Returns:
        Unit: libSBML unit object
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
        - model_specs (tuple): Level & Version of SBML model
        - identifier (str):    Identifier for the defined unit
        - name (str):          Full name of the defined unit
        - units (list):        All units the defined unit consists of
         
    Returns:
        UnitDefinition: libSBML unit definition object
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
        - model (libModel): Model loaded with libSBML
         
    Returns:
        list: List of libSBML UnitDefinitions
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
        - contained_unit_defs (list): List of libSBML UnitDefinition objects
    """
    for unit_def in contained_unit_defs:
        logging.info(unit_def.toXMLNode())


def print_remaining_UnitDefinitions(model: libModel, list_of_fba_units: list[UnitDefinition]):
    """Prints UnitDefinitions from the model that were removed as these were not contained in the list_of_fba_units

    Args:
        - model (libModel): Model loaded with libSBML
        - list_of_fba_units (list):  List of libSBML UnitDefinitions  
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
def add_fba_units(model: libModel):
    """| Adds:
       |     - mmol per gDW per h
       |     - mmol per gDW 
       |     - hour (h)
       |     - femto litre (fL)
       |
       | to the list of unit definitions (needed for FBA)

    Args:
        - model (libModel): Model loaded with libSBML
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
        - model (libModel): Model loaded with libSBML
    """ 
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


#--------------------------------------- Function to set units fro parameters -----------------------------------------#
def set_units(model: libModel):
    """Sets units of parameters in model

    Args:
        - model (libModel): Model loaded with libSBML
    """
    for param in model.getListOfParameters(): # needs to be added to list of unit definitions aswell
        if any(
            (unit_id := re.fullmatch('mmol_per_gDW_per_hr?', unit.getId(), re.IGNORECASE)) 
            for unit in model.getListOfUnitDefinitions()
            ):
            if not (param.isSetUnits() and param.getUnits() == unit_id.group(0)):
                param.setUnits(unit_id.group(0))
            

#-------------------------- Functions to add default settings for compartments & metabolites --------------------------#            
def add_compartment_structure_specs(model: libModel):
    """| Adds the required specifications for the compartment structure
       | if not set (size & spatial dimension)
        
    Args:
        - model (libModel): Model loaded with libSBML
    """ 
    for compartment in model.getListOfCompartments():
      
        if not compartment.isSetSize():
            compartment.setSize(float('NaN'))
         
        if not compartment.isSetSpatialDimensions():
            compartment.setSpatialDimensions(3)
         
        if any(
            (unit_id := re.fullmatch('fL', unit.getId(), re.IGNORECASE)) for unit in model.getListOfUnitDefinitions()
            ):
            if not (compartment.isSetUnits() and compartment.getUnits() == unit_id.group(0)):
                compartment.setUnits(unit_id.group(0))
         
         
def set_initial_amount(model: libModel):
    """Sets initial amount to all metabolites if not already set or if initial concentration is not set
    
    Args:
        - model (libModel): Model loaded with libSBML
    """
    for species in model.getListOfSpecies():
      
      if not (species.isSetInitialAmount() or species.isSetInitialConcentration()):
         species.setInitialAmount(float('NaN'))
         

#--------------------------------- Function to add URIs from the IDs for GeneProducts ---------------------------------# 
def cv_ncbiprotein(gene_list, email, protein_fasta: str, lab_strain: bool=False):
    """Adds NCBI Id to genes as annotation

    Args:
        - gene_list (list): libSBML ListOfGenes
        - email (str): User Email to access the Entrez database
        - protein_fasta (str): The path to the CarveMe protein.fasta input file
        - lab_strain (bool): Needs to be set to True if strain was self-annotated
                           and/or the locus tags in the CarveMe input file should be kept   
    """
    Entrez.email = email
                    
    id2locus_name = None  # Needs to be initialised, otherwise UnboundLocalError: local variable 'id2locus_name' referenced before assignment          
    if protein_fasta:
        if protein_fasta.strip() != '': 
            id2locus_name = parse_fasta_headers(protein_fasta)
            id2locus_name.set_index('protein_id')
    
    genes_missing_annotation = []

    print('Setting CVTerms and removing notes for all genes:')
    for gene in tqdm(gene_list):
        if not gene.isSetMetaId():
            gene.setMetaId('meta_' + gene.getId())
        
        if (gene.getId()[2] == 'W'): #addition to work with KC-Na-01
            entry = gene.getId()
            entry = entry[2:-7] if '__46__' in entry else entry[2:-2] # Required for VMH models
            add_cv_term_genes(entry, 'REFSEQ', gene)
            add_cv_term_genes(entry, 'NCBI', gene, lab_strain)
            name, locus = search_ncbi_for_gpr(entry)
            gene.setName(name)
            gene.setLabel(locus)
        
        elif (gene.getId() != 'G_spontaneous') and (gene.getId() != 'G_Unknown'): # Has to be omitted as no additional data can be retrieved neither from NCBI nor the CarveMe input file
            if 'prot_' in gene.getId():
                id_string = gene.getId().split('prot_')[1].split('_')  # All NCBI CDS protein FASTA files have the NCBI protein identifier after 'prot_' in the FASTA identifier
                ncbi_id = id_string[0]  # If identifier contains no '_', this is full identifier
            else:
                id_string = gene.getId().removeprefix('G_').split('_')
                if 'peg' in id_string: 
                    genes_missing_annotation.append('_'.join(id_string))
                    continue
              
            if len(id_string) == 2: # Can be the case if ID is locus tag, for example
                genes_missing_annotation.append('_'.join(id_string))
                continue # Ignore locus tags as no valid identifiers
            if (len(id_string) > 2):  # Identifier contains '_'
            # Check that the second entry consists of a sequence of numbers -> Valid RefSeq identifier! 
            # (Needs to be changed if there are other gene idenitfiers used that could contain '_' & need to be handled differently)
                if re.fullmatch('^\d+\d+$', id_string[1], re.IGNORECASE):
                    ncbi_id = '_'.join(id_string[:2])  # Merge the first two parts with '_' as this is complete identifier
           
            # If identifier matches RefSeq ID pattern   
            if re.fullmatch('^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', ncbi_id, re.IGNORECASE):
                add_cv_term_genes(ncbi_id, 'REFSEQ', gene, lab_strain)
                add_cv_term_genes(ncbi_id, 'NCBI', gene, lab_strain)
                name, locus = search_ncbi_for_gpr(ncbi_id)

            # If identifier only contains numbers 
            # -> Get the corresponding data from the CarveMe input file
            elif re.fullmatch('^\d+$', ncbi_id, re.IGNORECASE):
                if id2locus_name is not None:
                    name, locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['name', 'locus_tag']].values[0]
                else: 
                    genes_missing_annotation.append(ncbi_id)
        
            # If identifier matches ncbiprotein ID pattern
            elif re.fullmatch('^(\w+\d+(\.\d+)?)|(NP_\d+)$', ncbi_id, re.IGNORECASE):
                add_cv_term_genes(ncbi_id, 'NCBI', gene, lab_strain)
                name, locus = search_ncbi_for_gpr(ncbi_id)
            
            # Catch all remaining cases that have no valid ID   
            else: 
                genes_missing_annotation.append(ncbi_id)
        
            # For lab strains use the locus tag from the annotation file   
            if lab_strain and id2locus_name is not None:
                locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['locus_tag']].values[0][0]
        
            if ncbi_id not in genes_missing_annotation:
                gene.setName(name)
                gene.setLabel(locus)
            
        gene.unsetNotes()
    if genes_missing_annotation:    
        logging.warning(f'The following {len(genes_missing_annotation)} genes have no annotation, name & label (locus tag): {genes_missing_annotation}')


#------------------- Functions to change the CURIE pattern/CVTerm qualifier & qualifier type --------------------------#
def get_set_of_curies(uri_list: list[str]) -> tuple[SortedDict[str: SortedSet[str]], list[str]]:
    """| Gets a list of URIs
       | & maps the database prefixes to their respective identifier sets
        
    Args:
        - uri_list (list[str]): List containing CURIEs
            
    Returns:
        tuple: Two dictionaries (1) & (2)
            (1) SortedDict: Sorted dictionary mapping database prefixes from the provided CURIEs to their respective identifier sets also provided by the CURIEs
            (2) list: List of CURIEs that are invalid according to bioregistry
    """
    curie_dict = SortedDict()
    prefix, identifier = None, None
    invalid_curies = []
    
    for uri in uri_list:
        
        # Extracts the CURIE part from the URI/IRI
        if MIRIAM in uri:
            extracted_curie = uri.split(MIRIAM)[1]
        else:
            extracted_curie = uri.split(OLD_MIRIAM)[1]

        curie = manager.parse_curie(extracted_curie) # Contains valid db prefix to identifiers pairs
        curie = list(curie) # Turn tuple into list to allow item assignment
        
        if curie[0]: # Prefix is valid but to have same result for same databases need to do a bit of own parsing
            if re.fullmatch('^biocyc$', curie[0], re.IGNORECASE):  # Check for biocyc to also add metacyc if possible
                # Always add META if BioCyc sub-datbase prefixes are missing
                curie = curie if curie[1].split(':')[0] in BIOCYC_TIER1_DATABASES_PREFIXES else [curie[0], f'META:{curie[1]}']

                if 'META' in curie[1]: 
                    if is_valid_identifier(*curie): # Get the valid BioCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f'{curie[0]}:{curie[1]}')

                    # Add the MetaCyc identifier additionally
                    curie[1] = curie[1].split('META:')[1] # Metacyc identifier comes after 'META:' in biocyc identifier
                    if re.search('^rxn-|-rxn$', curie[1], re.IGNORECASE):
                        curie[0] = 'metacyc.reaction'
                    else:
                        curie[0] = 'metacyc.compound'
            elif 'metacyc.' in curie[0]:
                if is_valid_identifier(*curie): # Get the valid MetaCyc identifier & Add to dictionary
                    prefix, identifier = normalize_parsed_curie(*curie)

                    if not curie_dict or (prefix not in curie_dict):
                        curie_dict[prefix] = SortedSet()
                    curie_dict[prefix].add(identifier)
                else:
                    invalid_curies.append(f'{curie[0]}:{curie[1]}')

                # Add the BioCyc identifier additionally
                curie = ['biocyc', f'META:{curie[1]}'] # Metacyc identifier comes after 'META:' in biocyc identifier
            elif re.fullmatch('^brenda$', curie[0], re.IGNORECASE): # Brenda & EC code is the same
                curie[0] = 'eccode'
        
        elif not curie[0]: # Need to do own parsing if prefix is not valid
            # Get CURIEs irrespective of pattern
            if '/' in extracted_curie:
                extracted_curie = extracted_curie.split('/')

                # Check for NaN identifiers
                if re.fullmatch('^nan$', extracted_curie[0], re.IGNORECASE) or re.fullmatch('^nan$', extracted_curie[1], re.IGNORECASE):
                    # Only return strings where the database prefix is 'NaN' but a possible identifier could be contained
                    if re.fullmatch('^nan$', extracted_curie[0], re.IGNORECASE) and not re.fullmatch('^nan$', extracted_curie[1], re.IGNORECASE): 
                        invalid_curies.append(f'{extracted_curie[0]}:{extracted_curie[1]}')
                    continue
                # Check for certain special cases
                if re.search('inchi', extracted_curie[0], re.IGNORECASE):  # Check for inchi as splitting by '/' splits too much
                    if re.fullmatch('^inchi$', extracted_curie[0], re.IGNORECASE):
                        curie = (extracted_curie[0].lower(), '/'.join(extracted_curie[1:len(extracted_curie)]))
                    elif re.fullmatch('^inchikey$', extracted_curie[0], re.IGNORECASE):
                        curie = (extracted_curie[0].lower(), extracted_curie[1])
                    else:
                        wrong_prefix = extracted_curie[0].split(':')
                        curie = (wrong_prefix[0], f'{wrong_prefix[1]}/{"/".join(extracted_curie[1:len(extracted_curie)])}')
                elif re.fullmatch('^brenda$', extracted_curie[0], re.IGNORECASE): # Brenda & EC code is the same
                    curie = ('eccode', extracted_curie[1])
                elif re.fullmatch('^biocyc$', extracted_curie[0], re.IGNORECASE):  # Check for biocyc to also add metacyc if possible
                    # Always add META if BioCyc sub-datbase prefixes are missing
                    extracted_curie[1] = extracted_curie[1] if extracted_curie[1].split(':')[0] in BIOCYC_TIER1_DATABASES_PREFIXES else f'META:{extracted_curie[1]}'
                    curie = ['biocyc', extracted_curie[1]]

                    if 'META' in curie[1]: 
                        if is_valid_identifier(*curie): # Get the valid BioCyc identifier & Add to dictionary
                            prefix, identifier = normalize_parsed_curie(*curie)

                            if not curie_dict or (prefix not in curie_dict):
                                curie_dict[prefix] = SortedSet()
                            curie_dict[prefix].add(identifier)
                        else:
                            invalid_curies.append(f'{curie[0]}:{curie[1]}')

                        # Add additionallly the MetaCyc identifier 
                        curie[1] = curie[1].split('META:')[1] # Metacyc identifier comes after 'META:' in biocyc identifier
                        if re.search('^rxn-|-rxn$', curie[1], re.IGNORECASE):
                            curie[0] = 'metacyc.reaction'
                        else:
                            curie[0] = 'metacyc.compound'
                elif 'metacyc.' in extracted_curie[0]:
                    curie = extracted_curie
                    if is_valid_identifier(*curie): # Get the valid MetaCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f'{curie[0]}:{curie[1]}')

                    # Add BioCyc identfier additionally
                    curie = ['biocyc', f'META:{curie[1]}'] # Metacyc identifier comes after 'META:' in biocyc identifier
                elif re.fullmatch('^chebi$', extracted_curie[0], re.IGNORECASE):
                    new_curie = extracted_curie[1].split(':')
                    curie = (new_curie[0].lower(), new_curie[1])
                elif re.search('^sbo:', extracted_curie[1], re.IGNORECASE): # Checks for old pattern of SBO term URIs ('MIRIAM/sbo/SBO:identifier')
                    curie = [extracted_curie[0], extracted_curie[1].split(':')[1]]
                else:
                    if re.fullmatch('^brenda$', extracted_curie[0], re.IGNORECASE) or re.fullmatch('^ec-code$', extracted_curie[0], re.IGNORECASE): # Brenda equals EC code, EC code in URI = ec-code
                        curie[0] = 'eccode'
                    else:
                        curie[0] = extracted_curie[0]

                    curie[1] = extracted_curie[1]

            elif ':' in extracted_curie:
                extracted_curie = extracted_curie.split(':')
                
                # Check for NaN identifiers
                if re.fullmatch('^nan$', extracted_curie[0], re.IGNORECASE) or re.fullmatch('^nan$', extracted_curie[1], re.IGNORECASE):
                    # Only return strings where the database prefix is 'NaN' but a possible identifier could be contained
                    if re.fullmatch('^nan$', extracted_curie[0], re.IGNORECASE) and not re.fullmatch('^nan$', extracted_curie[1], re.IGNORECASE): 
                        invalid_curies.append(f'{extracted_curie[0]}:{extracted_curie[1]}')
                    continue
                elif re.fullmatch('^biocyc$', extracted_curie[0], re.IGNORECASE):  # Check for biocyc to also add metacyc if possible
                    # Always add META if BioCyc sub-datbase prefixes are missing
                    extracted_curie[1] = extracted_curie[1] if extracted_curie[1].split(':')[0] in BIOCYC_TIER1_DATABASES_PREFIXES else f'META:{extracted_curie[1]}'
                    curie = ['biocyc', extracted_curie[1]]

                    if 'META' in curie[1]:
                        if is_valid_identifier(*curie): # Get the valid BioCyc identifier & Add to dictionary
                            prefix, identifier = normalize_parsed_curie(*curie)

                            if not curie_dict or (prefix not in curie_dict):
                                curie_dict[prefix] = SortedSet()
                            curie_dict[prefix].add(identifier)
                        else:
                            invalid_curies.append(f'{curie[0]}:{curie[1]}')

                        # Add MetaCyc identifier additionally
                        curie[1] = curie[1].split('META:')[1] # Metacyc identifier comes after 'META:' in biocyc identifier
                        if re.search('^rxn-|-rxn$', curie[1], re.IGNORECASE):
                            curie[0] = 'metacyc.reaction'
                        else:
                            curie[0] = 'metacyc.compound'
                elif 'metacyc.' in extracted_curie[0]:
                    curie = extracted_curie
                    if is_valid_identifier(*curie): # Get the valid MetaCyc identifier & Add to dictionary
                        prefix, identifier = normalize_parsed_curie(*curie)

                        if not curie_dict or (prefix not in curie_dict):
                            curie_dict[prefix] = SortedSet()
                        curie_dict[prefix].add(identifier)
                    else:
                        invalid_curies.append(f'{curie[0]}:{curie[1]}')

                    # Add BioCyc identifier additionally
                    curie = ['biocyc', f'META:{curie[1]}'] # Metacyc identifier comes after 'META:' in biocyc identifier
                else:
                    if re.fullmatch('^brenda$', extracted_curie[0], re.IGNORECASE) or re.fullmatch('^ec-code$', extracted_curie[0], re.IGNORECASE): # Brenda equals EC code, EC code in URI = ec-code
                        curie[0] = 'eccode'
                    else:
                        curie[0] = extracted_curie[0]

                    if re.fullmatch('^kegg.genes$', extracted_curie[0], re.IGNORECASE):
                        curie[1] = ':'.join(extracted_curie[1:len(extracted_curie)])
                    else:
                        curie[1] = extracted_curie[1]

        if is_valid_identifier(*curie): # Get all valid identifiers
            prefix, identifier = normalize_parsed_curie(*curie)
        else:
            
            if curie[0] == 'eccode':
                correct_id = curie[1] # EC number needs to have 4 places if splitted at the dots
                while len(correct_id.split('.')) < 4:
                    correct_id = f'{correct_id}.-'
                prefix, identifier = normalize_parsed_curie(curie[0], correct_id)
                # Add too long EC codes back in model BUT report as invalid CURIEs!
                if (len(correct_id.split('.')) > 4): invalid_curies.append(f'{prefix}:{identifier}')
            # Rhea identifier should only contain 5 numbers but added by CarveMe the Rhea identifier contains '#1'
            elif (curie[0] == 'rhea') and ('#' in curie[1]):
                prefix, identifier = normalize_parsed_curie(curie[0], curie[1].split('#')[0])
            else:
                invalid_curies.append(f'{curie[0]}:{curie[1]}')
                
        if prefix and identifier: # Check that a prefix & identifier pair was found!
            # Use prefix as key & the corresponding set of identifiers as values   
            if not curie_dict or (prefix not in curie_dict):
                curie_dict[prefix] = SortedSet()
            curie_dict[prefix].add(identifier)
            
    return curie_dict, invalid_curies

def generate_uri_set_with_specific_pattern(prefix2id: SortedDict[str: SortedSet[str]], new_pattern: bool) -> SortedSet[str]: 
    """Generate a set of complete URIs from the provided prefix to identifier mapping
        
    Args:
        - prefix2id (SortedDict[str: SortedSet[str]]): Dictionary containing a mapping from database prefixes to their respective identifier sets 
        - new_pattern (bool):                          True if new pattern is wanted, otherwise False
            
    Returns:
        SortedSet: Sorted set containing complete URIs
    """
    uri_set = SortedSet()
    
    if new_pattern:
        SEPARATOR = ':'
    else:
        SEPARATOR = '/'
    
    for prefix in prefix2id:
        current_prefix = prefix   
        
        for identifier in prefix2id.get(prefix):
            separator = SEPARATOR

            if re.search('o$', prefix, re.IGNORECASE):  # Ontologies seem only to work with new pattern!
                separator = ':'
                prefix = prefix.upper()
            
            elif re.fullmatch('^chebi$', current_prefix, re.IGNORECASE) and not new_pattern:  # The old pattern for chebi is different: Just adding '/' das NOT work!
                prefix = f'chebi/{current_prefix}'
                separator = ':'
                
            elif re.fullmatch('^biocyc$', prefix, re.IGNORECASE):  # Get identifier for biocyc
                prefix = f'biocyc{SEPARATOR}META'
                separator = ':'

            
            uri = MIRIAM + prefix + separator + identifier
            uri_set.add(uri)
            
    return uri_set


def generate_miriam_compliant_uri_set(prefix2id: SortedDict[str: SortedSet[str]]) -> SortedSet[str]: 
    """Generate a set of complete MIRIAM compliant URIs from the provided prefix to identifier mapping
        
    Args:
        - prefix2id (SortedDict[str: SortedSet[str]]): Dictionary containing a mapping from database prefixes to their respective identifier sets
            
    Returns:
        SortedSet: Sorted set containing complete URIs
    """
    uri_set = SortedSet()
    
    for prefix in prefix2id:    
        for identifier in prefix2id.get(prefix):
            uri = get_identifiers_org_iri(prefix, identifier)
            uri_set.add(uri)
            
    return uri_set


def add_uri_set(entity: SBase, qt, b_m_qt, uri_set: SortedSet[str]) -> list[str]:
    """Add a complete URI set to the provided CVTerm
        
    Args:
        - entity (SBase):               A libSBML SBase object like model, GeneProduct, etc.
        - qt:                           A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - b_m_qt:                       A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - uri_set (SortedSet[str]):     SortedSet containing URIs
    """        
    new_cvterm = generate_cvterm(qt, b_m_qt)
        
    for uri in uri_set:
        new_cvterm.addResource(uri)
            
    entity.addCVTerm(new_cvterm)

def improve_uri_per_entity(entity: SBase, bioregistry: bool, new_pattern: bool) -> tuple[list[str], list[str]]:
    """Helper function: Removes duplicates & changes pattern according to new_pattern

    Args:
        - entity (SBase):       A libSBML SBase object, either a model or an entity
        - bioregistry (bool):   Specifies whether the URIs should be changed with the help of bioregistry to be MIRIAM compliant or changed according to new or old pattern
        - new_pattern (bool):   True if new pattern is wanted, otherwise False
        
    Returns:
        tuple: Two lists (1) & (2)
            (1) list: List of all collected invalid annotations of one entity
            (2) list: List of all collected invalid CURIEs of one entity
    """
    not_miriam_compliant = []
    collected_invalid_curies = []
    pattern = f'{MIRIAM}|{OLD_MIRIAM}'
    cvterms = entity.getCVTerms()
    
    for cvterm in cvterms:
        tmp_list = []
        
        # Retrieve QualifierType & Biological/ModelQualifierType before resource is removed!
        current_qt = cvterm.getQualifierType()
                
        if current_qt == BIOLOGICAL_QUALIFIER:
            current_b_m_qt = cvterm.getBiologicalQualifierType()
        elif current_qt == MODEL_QUALIFIER:
            current_b_m_qt = cvterm.getModelQualifierType()
            
        current_uris = [cvterm.getResourceURI(i) for i in range(cvterm.getNumResources())]
    
        for cu in current_uris:
            if re.match(pattern, cu, re.IGNORECASE):  # If model contains identifiers without MIRIAM/OLD_MIRIAM these are kept 
                tmp_list.append(cu)
                # Remove all valid URIs to add these back later again
                # In case of EC numbers, remove all EC number URIs as EC numbers not containing four numbers are sometimes invalid and sometimes not
                if is_valid_identifier(*manager.parse_uri(cu)) or (manager.parse_uri(cu)[0] == 'eccode'): cvterm.removeResource(cu)
            else:
                not_miriam_compliant.append(cu)
            
        prefix2id, invalid_curies = get_set_of_curies(tmp_list)
        collected_invalid_curies.extend(invalid_curies)
        if bioregistry: uri_set = generate_miriam_compliant_uri_set(prefix2id)
        else: uri_set = generate_uri_set_with_specific_pattern(prefix2id, new_pattern)
        add_uri_set(entity, current_qt, current_b_m_qt, uri_set)
    
    return not_miriam_compliant, collected_invalid_curies

def improve_uris(entities: SBase, bioregistry: bool, new_pattern: bool) -> tuple[dict[str:list[str]], dict[str:list[str]]]:
    """Removes duplicates & changes pattern according to bioregistry or new_pattern
    
    Args:
        - entities (SBase):     A libSBML SBase object, either a model or a list of entities
        - bioregistry (bool):   Specifies whether the URIs should be changed with the help of bioregistry to be MIRIAM compliant or changed according to new or old pattern
        - new_pattern (bool):   True if new pattern is wanted, otherwise False

    Returns:
        tuple: Two dictionnaries (1) & (2)
            (1) dictionary: Mapping of entity identifier to list of corresponding not MIRIAM compliant annotations 
            (2) dictionary: Mapping of entity identifier to list of corresponding invalid CURIEs
    """
    entity2not_miriam = {}
    entity2invalid_curies = {}

    if type(entities) == libModel:  # Model needs to be handled like entity!
        not_miriam_compliant, invalid_curies = improve_uri_per_entity(entities, bioregistry, new_pattern)
        if not_miriam_compliant: entity2not_miriam[entities.getId()] = not_miriam_compliant
        if invalid_curies: entity2invalid_curies[entities.getId()] = invalid_curies
    
    else: 
        for entity in tqdm(entities):
            not_miriam_compliant, invalid_curies = improve_uri_per_entity(entity, bioregistry, new_pattern)
            if not_miriam_compliant: entity2not_miriam[entity.getId()] = not_miriam_compliant
            if invalid_curies: entity2invalid_curies[entity.getId()] = invalid_curies
            
            if type(entity) == UnitDefinition:
                for unit in entity.getListOfUnits():  # Unit needs to be handled within ListOfUnitDefinition
                    not_miriam_compliant, invalid_curies = improve_uri_per_entity(unit, bioregistry, new_pattern)
                    if not_miriam_compliant: entity2not_miriam[unit.getId()] = not_miriam_compliant
                    if invalid_curies: entity2invalid_curies[unit.getId()] = invalid_curies
                    
    return entity2not_miriam, entity2invalid_curies


def polish_annotations(model: libModel, bioregistry: bool, new_pattern: bool, filename: str) -> libModel:
    """| Polishes all annotations in a model such that no duplicates are present 
       | & the same pattern is used for all CURIEs
        
    Args:
        - model (libModel):     Model loaded with libSBML
        - bioregistry (bool):   Specifies whether the URIs should be changed with the help of bioregistry to be MIRIAM compliant or changed according to new or old pattern
        - new_pattern (bool):   True if new pattern is wanted, otherwise False
        - filename (str):       Path to output file for invalid CURIEs detected by improve_uris
        
    Returns:
        libModel: libSBML model with polished annotations
    """
    list_of_entity2not_miriam, list_of_entity2invalid_curies = [], []
    listOf_dict = {
        'model': model,
        'compartment': model.getListOfCompartments(),
        'metabolite': model.getListOfSpecies(),
        'parameter': model.getListOfParameters(),
        'reaction': model.getListOfReactions(),
        'unit definition': model.getListOfUnitDefinitions(),
        }

    if model.isPackageEnabled('fbc'):
        listOf_dict['gene product'] = model.getPlugin('fbc').getListOfGeneProducts()
    
    if model.isPackageEnabled('groups'):
        listOf_dict['group'] = model.getPlugin('groups').getListOfGroups()

    # Adjust annotations in model
    for listOf in listOf_dict:
        print(f'Polish {listOf} annotations...')
        entity2not_miriam, entity2invalid_curies = improve_uris(listOf_dict[listOf], bioregistry, new_pattern)
        list_of_entity2not_miriam.append(entity2not_miriam)
        list_of_entity2invalid_curies.append(entity2invalid_curies)
        
    all_entity2not_miriam = reduce(lambda d1, d2: {**d1, **d2}, list_of_entity2not_miriam)
    all_entity2invalid_curies = reduce(lambda d1, d2: {**d1, **d2}, list_of_entity2invalid_curies)
        
    if all_entity2not_miriam:
        miriam_filename = f'{filename}_invalid_annotations_{str(date.today().strftime("%Y%m%d"))}.tsv'
        logging.warning(f'In the provided model {model.getId()} for {len(all_entity2not_miriam)} entities invalid annotations were detected. ' +
                     f'These invalid annotations are saved to {miriam_filename}')
        not_miriam_compliant_df = parse_dict_to_dataframe(all_entity2not_miriam)
        not_miriam_compliant_df.columns = ['entity', 'not_miriam']
        not_miriam_compliant_df.to_csv(miriam_filename, sep='\t')
    
    if all_entity2invalid_curies: 
        curies_filename = f'{filename}_invalid_curies_{str(date.today().strftime("%Y%m%d"))}.tsv'      
        logging.warning(f'In the provided model {model.getId()} for {len(all_entity2invalid_curies)} entities invalid CURIEs were detected. ' +
                     f'These invalid CURIEs are saved to {curies_filename}')      
        invalid_curies_df = parse_dict_to_dataframe(all_entity2invalid_curies)
        invalid_curies_df.columns = ['entity', 'invalid_curie']
        invalid_curies_df[['prefix', 'identifier']] = invalid_curies_df.invalid_curie.str.split(':', n=1, expand = True) # Required for identifiers that aso contain a ':'
        invalid_curies_df = invalid_curies_df.drop('invalid_curie', axis=1)
        invalid_curies_df.to_csv(curies_filename, sep='\t')
    
    return model


def change_qualifier_per_entity(entity: SBase, new_qt, new_b_m_qt, specific_db_prefix: str=None) -> list:
    """Updates Qualifiers to be MIRIAM compliant for an entity

    Args:
        - entity (SBase): A libSBML SBase object like model, GeneProduct, etc.
        - new_qt (Qualifier): A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - new_b_m_qt (QualifierType): A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - specific_db_prefix (str): Has to be set if only for a specific database the qualifier type should be changed. Can be 'kegg.genes', 'biocyc', etc.

    Returns:
        list: CURIEs that are not MIRIAM compliant
    """
    not_miriam_compliant = []
    pattern = f'{MIRIAM}|{OLD_MIRIAM}'
    cvterms = entity.getCVTerms()

    #for i in range(len(cvterms)):
    for cvterm in cvterms:
        tmp_set = SortedSet()
        sbo_set = SortedSet()
        #cvterm = cvterms.get(i)
        
        # include check for reaction and unit definition
        # if entity == Reaction or entity == UnitDefinition:
        # print(cvterm.getBiologicalQualifierType())
        if cvterm.getBiologicalQualifierType() == 9:  # 9 = BQB_OCCURS_IN (Reaction), Check for reactions with occursIn
            logging.info(f'CVTerm for {Fore.LIGHTYELLOW_EX}{str(entity)}{Style.RESET_ALL}' +
                  f' is left as {Fore.LIGHTYELLOW_EX}{BiolQualifierType_toString(cvterm.getBiologicalQualifierType())}{Style.RESET_ALL}')
        
        elif cvterm.getModelQualifierType() == 1:  # 1 = BQM_IS_DESCRIBED_BY (UnitDefinition), Check for UnitDefinitions with isDescribedBy
            logging.info(f'CVTerm for {Fore.LIGHTYELLOW_EX}{str(entity)}{Style.RESET_ALL}' + 
                  f' is left as {Fore.LIGHTYELLOW_EX}{ModelQualifierType_toString(cvterm.getModelQualifierType())}{Style.RESET_ALL}')
        
        else:
            current_curies = [cvterm.getResourceURI(j) for j in range(cvterm.getNumResources())]
        
            for cc in current_curies:
                    
                current_curie = None
                    
                if (specific_db_prefix != None) and (specific_db_prefix != ''):
                    if specific_db_prefix in cc:
                        current_curie = cc
                else:
                    current_curie = cc
                    
                if (current_curie) and re.match(pattern, current_curie, re.IGNORECASE):  # If model contains identifiers without MIRIAM/OLD_MIRIAM these are kept
                    if re.search('sbo:', current_curie, re.IGNORECASE): sbo_set.add(current_curie)
                    else: tmp_set.add(current_curie)
                    cvterm.removeResource(current_curie)
                else:
                    not_miriam_compliant.append(current_curie)
            
            if sbo_set: add_uri_set(entity, BIOLOGICAL_QUALIFIER, BQB_HAS_PROPERTY, sbo_set)
            add_uri_set(entity, new_qt, new_b_m_qt, tmp_set)
            #cvterms.remove(i)
                
    if not_miriam_compliant:
        return not_miriam_compliant


def change_qualifiers(model: libModel, entity_type: str, new_qt, new_b_m_qt, specific_db_prefix: str = None) -> libModel:
    """Updates Qualifiers to be MIRIAM compliant for an entity type of a given model 

    Args:
        - model (libModel):   Model loaded with libSBML
        - entity_type (str): Any string of the following: model|compartment|metabolite|parameter|reaction|unit definition|unit|gene product|group
        - new_qt (Qualifier): A libSBML qualifier type: BIOLOGICAL_QUALIFIER|MODEL_QUALIFIER
        - new_b_m_qt (QualifierType): A libSBML biological or model qualifier type like BQB_IS|BQM_IS
        - specific_db_prefix (str): Has to be set if only for a specific database the qualifier type should be changed. Can be 'kegg.genes', 'biocyc', etc.
    
    Returns:
        libModel: Model with changed qualifier for given entity type
    """
    not_miriam_compliant = []
    listOf_dict = {
        'model': model,
        'compartment': model.getListOfCompartments(),
        'metabolite': model.getListOfSpecies(),
        'parameter': model.getListOfParameters(),
        'reaction': model.getListOfReactions(),
        'unit definition': model.getListOfUnitDefinitions(),
        }
    
    if model.isPackageEnabled('fbc'):
        listOf_dict['gene product'] = model.getPlugin('fbc').getListOfGeneProducts()
    
    if model.isPackageEnabled('groups'):
        listOf_dict['group'] = model.getPlugin('groups').getListOfGroups()
        
    if entity_type == 'model':  # Model needs to be handled like entity!
        not_miriam_compliant = change_qualifier_per_entity(listOf_dict.get('model'), new_qt, new_b_m_qt, specific_db_prefix)
        
    elif entity_type == 'unit':
        for unit in listOf_dict.get('unit definition'):  # Unit needs to be handled within ListOfUnitDefinition
            not_miriam_compliant = change_qualifier_per_entity(unit, new_qt, new_b_m_qt, specific_db_prefix)
        
    else:
        try: 
            for entity in tqdm(listOf_dict.get(entity_type)):
                not_miriam_compliant = change_qualifier_per_entity(entity, new_qt, new_b_m_qt, specific_db_prefix)
        except(TypeError):
            logging.info('The entity ' +  entity_type + ' is not present in ' + model.getId())        
        
    if not_miriam_compliant:         
        logging.warning(f'The following {len(not_miriam_compliant)} entities are not MIRIAM compliant: {not_miriam_compliant}')
    
    return model


def change_all_qualifiers(model: libModel, lab_strain: bool) -> libModel:
    """Wrapper function to change qualifiers of all entities at once

    Args:
        - model (libModel): Model loaded with libSBML
        - lab_strain (bool): True if the strain was sequenced in a local lab

    Returns:
        libModel: Model with all qualifiers updated to be MIRIAM compliant
    """
    
    entity_list_mod = ['model',
                    'unit definition',
                    'unit']
    for entity in entity_list_mod:
        print(f'Change {str(entity)} qualifiers...')
        model = change_qualifiers(model, entity, MODEL_QUALIFIER, BQM_IS)
    
    entity_list = ['compartment',
                   'metabolite',
                   'parameter',
                   'reaction',
                   'gene product',
                   'group']
    for entity in entity_list:
        print(f'Change {str(entity)} qualifiers...')
        if lab_strain and entity == 'gene product':
            model = change_qualifiers(model, 'gene product', BIOLOGICAL_QUALIFIER, BQB_IS_HOMOLOG_TO)
        else:
            model = change_qualifiers(model, entity, BIOLOGICAL_QUALIFIER, BQB_IS)
        
    return model


#--------------------------------------------------- Main function ----------------------------------------------------#
def polish(model: libModel, email: str, id_db: str, protein_fasta: str, lab_strain: bool, path: str) -> libModel: 
    """| Completes all steps to polish a model
       | (Tested for models having either BiGG or VMH identifiers.)

    Args:
        - model (libModel): model loaded with libSBML
        - email (str): E-mail for Entrez
        - id_db (str): Main database where identifiers in model come from
        - protein_fasta (str): File used as input for CarveMe
        - lab_strain (bool): True if the strain was sequenced in a local lab
        - path (str): Output path for incorrect annotations file(s)
    
    Returns:
        libModel: Polished libSBML model
    """
    colorama_init(autoreset=True)
    
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
    
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    gene_list = model.getPlugin('fbc').getListOfGeneProducts()

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
    cv_ncbiprotein(gene_list, email, protein_fasta, lab_strain)
    
    ### set boundaries and constant ###
    polish_entities(metab_list, metabolite=True)
    polish_entities(reac_list, metabolite=False)

    
    ### MIRIAM compliance of CVTerms ###
    print('Remove duplicates & transform all CURIEs to the new identifiers.org pattern (: between db and ID):')
    filename = f'{path}{model.getId()}'
    polish_annotations(model, True, True, filename)
    print('Changing all qualifiers to be MIRIAM compliant:')
    change_all_qualifiers(model, lab_strain)
    
    return model