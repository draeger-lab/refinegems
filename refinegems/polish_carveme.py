#!/usr/bin/env python
""" Can be used to polish a model created with CarveMe v.1.5.1

The newer version of CarveMe leads to some irritations in the model, these scripts enable for example the addition of BiGG Ids to the annotations as well as a correct formatting of the annotations.
"""

from libsbml import *
from Bio import Entrez, SeqIO
from tqdm.auto import tqdm
from refinegems.cvterms import add_cv_term_metabolites, add_cv_term_reactions, add_cv_term_genes, metabol_db_dict, reaction_db_dict
from refinegems.load import write_to_file

def unset_ann(entity_list):
    """removes "empty" rdf bags from model

    Args:
        entity_list (list): libSBML ListOfSpecies or ListOfReactions
    """
    for entity in entity_list:
        ann = entity.getAnnotationString()
        if '<bqbiol:is xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">' in ann:
            entity.unsetAnnotation()
            
def add_bigg_metab(entity_list):
    """adds the BiGG Id of metabolites as URL to the annotation field

    Args:
        entity_list (list): libSBML ListOfSpecies
    """
    for entity in entity_list:
        bigg_id = entity.getId()
        bigg_id = bigg_id[2:] 
        bigg_id = bigg_id[:-2]
        
        add_cv_term_metabolites(bigg_id, 'BIGG', entity)
        # url = "https://identifiers.org/bigg.metabolite/" + bigg_id
        # cv = CVTerm()
        # cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        # cv.setBiologicalQualifierType(BQB_IS)
        # cv.addResource(url)
        # entity.addCVTerm(cv, False)
            
def add_bigg_reac(entity_list):
    """adds the BiGG Id of reactions as URL to the annotation field

    Args:
        entity_list (list): libSBML ListOfReactions
    """
    for entity in entity_list:
        bigg_id = entity.getId()
        if bigg_id != 'Growth':
            bigg_id = bigg_id[2:]
            add_cv_term_reactions(bigg_id, 'BIGG', entity)
            
        # url = "https://identifiers.org/bigg.reaction/" + bigg_id
        # cv = CVTerm()
        # cv.setQualifierType(BIOLOGICAL_QUALIFIER)
        # cv.setBiologicalQualifierType(BQB_IS)
        # cv.addResource(url)
        # entity.addCVTerm(cv, False)

def cv_notes_metab(species_list):
    """checks the notes field for information which should be in the annotation field
       removes entry from notes and adds it as URL to the CVTerms of a metabolite

    Args:
        species_list (list): libSBML ListOfSpecies
    """

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
                            add_cv_term_metabolites(entry.strip(), db, species)
                    else:
                        add_cv_term_metabolites(fill_in, db, species)
                        
        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)
        
        ####Adding new, shortened notes
        new_notes = ' '.join([str(elem)+'\n' for elem in notes_list])
        species.unsetNotes()
        species.setNotes(new_notes)
        #print(species.getNotesString())
                    
                    
def cv_notes_reac(reaction_list):
    """checks the notes field for information which should be in the annotation field
       removes entry from notes and adds it as URL to the CVTerms of a reaction

    Args:
        reaction_list (list): libSBML ListOfReactions
    """

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
                            add_cv_term_reactions(entry.strip(), db, reaction)
                    else:
                        add_cv_term_reactions(fill_in, db, reaction)
                        
        for elem in notes_string:
            if elem not in elem_used and elem not in notes_list:
                notes_list.append(elem)
        
        ####Adding new, shortened notes
        new_notes = ' '.join([str(elem)+'\n' for elem in notes_list])
        reaction.unsetNotes()
        reaction.setNotes(new_notes)
        
def polish_entities(entity_list, metabolite):
    """removes Notes field from entity
       sets boundary condition and constant if not set for a metabolite

    Args:
        entity_list (list): libSBML ListOfSpecies or ListOfReactions
        metabolite (boolean): flag to determine whether entity = metabolite
    """
    for entity in entity_list:
        entity.unsetNotes()
        if metabolite: # some polishing
            if not entity.getBoundaryCondition():
                entity.setBoundaryCondition(False)
            if not entity.getConstant():
                entity.setConstant(False)
                
def set_units(model):
    for param in model.getListOfParameters():
        if param.isSetUnits() == False:
            param.setUnits('mmol_per_gDW_per_hr')
                
def cv_ncbiprotein(gene_list, email):
    Entrez.email = email
    
    def get_name_locus_tag(ncbi_id):
        handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="gbwithparts", retmode='text')
        records = SeqIO.parse(handle, "gb")

        for i,record in enumerate(records):
            for feature in record.features:
                if feature.type == "CDS":
                    return record.description, feature.qualifiers["locus_tag"][0]

    print('Setting CVTerms and removing notes for all genes:')
    for gene in tqdm(gene_list):
        id_string = gene.getId().split('_')
        
        for entry in id_string:
            if (entry[:1] == 'Q'): #maybe change this to user input
                add_cv_term_genes(entry, 'NCBI', gene)
                name, locus = get_name_locus_tag(entry)
                gene.setName(name)
                gene.setLabel(locus)

        gene.unsetNotes()
            
def polish_carveme(model, new_filename, email):
    """completes all steps to polish a model created with CarveMe

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    metab_list = model.getListOfSpecies()
    reac_list = model.getListOfReactions()
    gene_list = model.getPlugin('fbc').getListOfGeneProducts()
    
    set_units(model)
    unset_ann(metab_list)
    unset_ann(reac_list)
    add_bigg_metab(metab_list)
    add_bigg_reac(reac_list)
    cv_notes_metab(metab_list)
    cv_notes_reac(reac_list)
    cv_ncbiprotein(gene_list, email)
    polish_entities(metab_list, metabolite=True)
    polish_entities(reac_list, metabolite=False)
    
    write_to_file(model, new_filename)