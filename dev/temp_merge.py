import re

from cobra.io.sbml import _f_specie, _f_reaction
from libsbml import Species, ListOfSpecies, ListOfReactions
from typing import Union

from refinegems.utility.cvterms import PREFIX2DB_METABS, PREFIX2DB_REACS, add_cv_term_metabolites, add_cv_term_reactions
from refinegems.utility.io import load_a_table_from_database
from refinegems.utility.util import DB2REGEX

def add_entity_temp(entity_list: Union[ListOfSpecies, ListOfReactions], id_db_prefix:str) -> None:

    # Set-up default variables
    try:
        db_pattern = DB2REGEX[id_db_prefix]
    except KeyError as e:
        print(f'{e}\n id_db_prefix must be one of the valid prefixes in https://bioregistry.io/.')
        return
        
    # Set-up case-dependent variables
    match entity_list:
        case ListOfSpecies():
            bigg_db = 'bigg_metabolites'
            bigg_id_type = 'universal_bigg_id'
            add_cv_term = add_cv_term_metabolites
            id_handler = _f_specie
            prefix2db = PREFIX2DB_METABS
        case ListOfReactions():
            bigg_db = 'bigg_reactions'
            bigg_id_type = 'id'
            add_cv_term = add_cv_term_reactions
            id_handler = _f_reaction
            prefix2db = PREFIX2DB_REACS
        case _:
            raise TypeError(f'Unsupported type for entity_list {type(entity_list)}. Must be ListOfSpecies or ListOfReactions.')

    # Get BiGG IDs for VMH ID == BiGG ID validation
    if 'vmh' in id_db_prefix.lower():
        bigg_ids = load_a_table_from_database(f'SELECT bigg_id FROM {bigg_db}')
        bigg_ids = set(bigg_ids[bigg_id_type].tolist())

    # Run operation
    for entity in entity_list:
        
        # Get ID
        current_id = entity.getId()

        # Use current_id as metaid if no metaid is present   
        if not entity.isSetMetaId():
            entity.setMetaId(f'meta_{current_id}')

        if re.fullmatch(db_pattern, current_id, re.IGNORECASE):

            # Remove prefix
            current_id = id_handler(current_id)
        
            # Unset annotations if no CV terms exist
            if entity.getNumCVTerms() == 0:
                entity.unsetAnnotation()

            # Get ID for annotation
            if isinstance(entity, Species): # Remove compartment suffix
                id_for_anno = re.sub(f'(_|\[){entity.getCompartment()}\]?$', '', current_id)
            else:
                id_for_anno = current_id

            # Add ID as URI to annotation
            add_cv_term(id_for_anno, prefix2db[id_db_prefix], entity)

            # Add BiGG ID to annotation, additionally, if VMH and valid BiGG ID
            if ('vmh' in id_db_prefix.lower()) and (id_for_anno in bigg_ids):
                    add_cv_term(id_for_anno, prefix2db['BIGG'], entity)
