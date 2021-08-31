#!/usr/bin/env python
""" Can be used to polish a model created with CarveMe v.1.5.1

The newer version of CarveMe leads to some irritations in the model, these scripts enable for example the addition of BiGG Ids to the annotations as well as a correct formatting of the annotations.
"""

from libsbml import *

def unset_ann(entity_list):
    for entity in entity_list:
        ann = entity.getAnnotationString()
        if '<bqbiol:is xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">' in ann:
            entity.unsetAnnotation()
            
def add_bigg_metab(entity_list):
    for entity in entity_list:
        bigg_id = entity.getId()
        bigg_id = bigg_id[2:] 
        bigg_id = bigg_id[:-2]
        url = "https://identifiers.org/bigg.metabolite/" + bigg_id
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(url)
        entity.addCVTerm(cv, False)
            
def add_bigg_reac(entity_list):
    for entity in entity_list:
        bigg_id = entity.getId()
        if bigg_id != 'Growth':
            bigg_id = bigg_id[2:] 
        url = "https://identifiers.org/bigg.reaction/" + bigg_id
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource(url)
        entity.addCVTerm(cv, False)

def cv_notes_metab(species_list):
    metabol_db_dict = {'BIGG': 'bigg.metabolite:', 'BRENDA': 'brenda:', 'CHEBI': 'chebi:','INCHI': 'inchi:', 'KEGG': 'kegg.compound:', 'METACYC': 'metacyc.compound:','MXNREF': 'metanetx.chemical:', 'SEED':'seed.compound:', 'UPA': 'unipathway.compound:', 'HMDB': 'hmdb:', 'REACTOME': 'reactome:'}
    
    def add_cv_term_from_notes_species(entry, db_id, metab):
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource('https://identifiers.org/'+ metabol_db_dict[db_id]+entry)
        metab.addCVTerm(cv)

    for species in species_list:
        notes_list=[]
        elem_used=[]
        notes_string = species.getNotesString().split('\n')
        for elem in notes_string:
            for db in metabol_db_dict.keys():
                if '<p>'+db in elem:
                    elem_used.append(elem)
                    #print(elem.strip()[:-4].split(': ')[1])
                    fill_in = elem.strip()[:-4].split(': ')[1]
                    if (';')  in fill_in and db != 'INCHI':
                        entries = fill_in.split(';')
                        for entry in entries:
                            add_cv_term_from_notes_species(entry.strip(), db, species)
                    else:
                        add_cv_term_from_notes_species(fill_in, db, species)
                        
        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)
        
        ####Adding new, shortened notes
        new_notes = ' '.join([str(elem)+'\n' for elem in notes_list])
        species.unsetNotes()
        species.setNotes(new_notes)
        #print(species.getNotesString())
                    
                    
def cv_notes_reac(reaction_list):
    reaction_db_dict = {'BIGG': 'bigg.reaction/', 'BRENDA': 'brenda/', 'KEGG': 'kegg.reaction/', 'METACYC': 'metacyc.reaction/', 'MetaNetX': 'metanetx.reaction/', 'SEED':'seed.reaction/','UPA': 'unipathway.reaction/', 'HMDB': 'hmdb/', 'REACTOME': 'reactome/', 'RHEA': 'rhea/','EC': 'ec-code/'}
    
    def add_cv_term_from_notes_reactions(entry, db_id, reac):
        cv = CVTerm()
        cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(BQB_IS)
        cv.addResource('https://identifiers.org/'+ reaction_db_dict[db_id]+entry)
        reac.addCVTerm(cv)

    for reaction in reaction_list:
        notes_list=[]
        elem_used=[]
        notes_string = reaction.getNotesString().split('\n')
        
        for elem in notes_string:
            for db in reaction_db_dict.keys():
                if '<p>'+db in elem:
                    elem_used.append(elem)
                    #print(elem.strip()[:-4].split(': ')[1])
                    fill_in = elem.strip()[:-4].split(': ')[1]
                    if (';')  in fill_in:
                        entries = fill_in.split(';')
                        for entry in entries:
                            add_cv_term_from_notes_reactions(entry.strip(), db, reaction)
                    else:
                        add_cv_term_from_notes_reactions(fill_in, db, reaction)
                        
        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)
        
        ####Adding new, shortened notes
        new_notes = ' '.join([str(elem)+'\n' for elem in notes_list])
        reaction.unsetNotes()
        reaction.setNotes(new_notes)
        
def polish_entities(entity_list, metabolite):
    for entity in entity_list:
        entity.unsetNotes()
        if metabolite: # some polishing
            if not entity.getBoundaryCondition():
                entity.setBoundaryCondition(False)
            if not entity.getConstant():
                entity.setConstant(False)
                
def write_to_file(model, new_filename):
    """Writes modified model to new file

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    new_document = model.getSBMLDocument()
    writeSBMLToFile(new_document, new_filename)
    print("Polished model written to " + new_filename)
            
def polish_carveme(model, new_filename):
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    
    unset_ann(metab_list)
    unset_ann(reac_list)
    add_bigg_metab(metab_list)
    add_bigg_reac(reac_list)
    cv_notes_metab(metab_list)
    cv_notes_reac(reac_list)
    polish_entities(metab_list, metabolite=True)
    polish_entities(reac_list, metabolite=False)
    
    write_to_file(model, new_filename)