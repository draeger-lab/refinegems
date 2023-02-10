#!/usr/bin/env python
""" Helper module to work with annotations (CVTerms)

Stores dictionaries which hold information the identifiers.org syntax, has functions to add CVTerms to different entities and parse CVTerms.
"""
from libsbml import BIOLOGICAL_QUALIFIER, BQB_IS, CVTerm, MODEL_QUALIFIER, BQM_IS, BQM_IS_DERIVED_FROM, BQM_IS_DESCRIBED_BY, Unit

metabol_db_dict = {
                   'BIGG': 'bigg.metabolite:',
                   'BIOCYC': 'metacyc.compound:',
                   'BioCyc': 'metacyc.compound:',
                   'BRENDA': 'brenda:',
                   'CHEBI': 'CHEBI:',
                   'ChEBI': 'CHEBI:',
                   'HMDB': 'hmdb:HMDB',
                   'Human Metabolome Database': 'hmdb:HMDB',
                   'INCHI': 'inchi:',
                   'InChI': 'inchi:',
                   'InChI-Key': 'inchikey:',
                   'KEGG': 'kegg.compound:',
                   #'KEGG Compound': 'kegg.compound:',
                   'METACYC': 'metacyc.compound:',
                   'MXNREF': 'metanetx.chemical:',
                   'MetaNetX': 'metanetx.chemical:',
                   'PUBCHEM': 'pubchem.compound:',
                   'REACTOME': 'reactome:',
                   'Reactome': 'reactome:',
                   'SEED': 'seed.compound:',
                   #'UPA': 'unipathway.compound:', #closed due to financial issues
                   #'UniPathway Compound': 'unipathway.compound:'
                   'VMH': 'vmhmetabolite:'
                  }

reaction_db_dict = {
                    'BIGG': 'bigg.reaction:',
                    'BioCyc': 'metacyc.reaction:',
                    'BRENDA': 'brenda:',
                    'EC': 'ec-code:',
                    'HMDB': 'hmdb:HMDB',
                    'KEGG': 'kegg.reaction:',
                    'METACYC': 'metacyc.reaction:',
                    'MXNREF': 'metanetx.reaction:',
                    'MetaNetX': 'metanetx.reaction:',
                    'REACTOME': 'reactome:',
                    'Reactome': 'reactome:',
                    'RHEA': 'rhea:',
                    'SEED': 'seed.reaction:',
                    #'UPA': 'unipathway.reaction:',
                    #'UniPathway Reaction': 'unipathway.reaction:'
                    'VMH': 'vmhreaction:'
                   }

gene_db_dict = {
                'KEGG': 'kegg.genes:',
                'NCBI': 'ncbiprotein:',
                'REFSEQ': 'refseq:',
                'UNIPROT': 'uniprot:'
               }

pathway_db_dict = {'KEGG': 'kegg.pathway:'}


def add_cv_term_units(unit_id: str, unit: Unit, relation: int):
    '''Adds CVTerm to a unit

       Params:
         unit_id (string):    ID to add as URI to annotation
         unit (libSBML Unit): Unit to add CVTerm to
         relation (str):      Provides model qualifier to be added
    '''
    cv = CVTerm()
    cv.setQualifierType(MODEL_QUALIFIER)
    cv.setModelQualifierType(relation)
    
    if relation == BQM_IS_DESCRIBED_BY:
      cv.addResource(f'https://identifiers.org/{unit_id}')
    else:
      cv.addResource(f'https://identifiers.org/UO:{unit_id}')
    unit.addCVTerm(cv)


def add_cv_term_metabolites(entry, db_id, metab):
    """Add CVTerm to a metabolite

    Args:
        entry (string): id to add as annotation
        db_id (string): database to which entry belongs
        metab (libsbml-species): metabolite to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + metabol_db_dict[db_id] + entry)
    metab.addCVTerm(cv)


def add_cv_term_reactions(entry, db_id, reac):
    """Add CVTerm to a reaction

    Args:
        entry (string): id to add as annotation
        db_id (string): database to which entry belongs
        reac (libsbml-reaction): reaction to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource(
        'https://identifiers.org/' +
        reaction_db_dict[db_id] +
        entry)
    reac.addCVTerm(cv)


def add_cv_term_genes(entry, db_id, gene):
    """Add CVTerm to a gene

    Args:
        entry (string): id to add as annotation
        db_id (string): database to which entry belongs
        gene (libsbml-gene): gene to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + gene_db_dict[db_id] + entry)
    gene.addCVTerm(cv)


def add_cv_term_pathways(entry, db_id, entity):
    """Add CVTerm to a groups pathway

    Args:
        entry (string): id to add as annotation
        db_id (string): database to which entry belongs
        entity (libsbml-group): pathway to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_IS)
    cv.addResource('https://identifiers.org/' + pathway_db_dict[db_id] + entry)
    entity.addCVTerm(cv)
    

def add_cv_term_pathways_to_entity(entry, db_id, entity):
    """Add CVTerm to a entity as OCCURS IN pathway

    Args:
        entry (string): id to add as annotation
        db_id (string): database to which entry belongs
        entity (libsbml-group): entity to add CVTerm to
    """
    cv = CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_OCCURS_IN)
    cv.addResource('https://identifiers.org/' + pathway_db_dict[db_id] + entry)
    entity.addCVTerm(cv)


def parse_id_from_cv_term(entity, db_id):
    """extract id for a specific database from CVTerm

    Args:
        entity (libsbml Object): Species, Reaction, Gene
        db_id (string): database of interest

    Returns:
        list: ids of entity belonging to db_id
    """
    num_cvs = entity.getNumCVTerms()
    all_ids = []
    for i in range(0, num_cvs):
        ann_string = entity.getCVTerm(i)
        num_res = ann_string.getNumResources()
        ids = [ann_string.getResourceURI(r).split(
            '/')[-1] for r in range(0, num_res) if str(db_id) in ann_string.getResourceURI(r)]
        ids = [id_string.split(':')[-1] for id_string in ids if ':' in id_string]
        all_ids.extend(ids)

    return all_ids
