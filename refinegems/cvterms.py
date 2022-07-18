#!/usr/bin/env python
""" Helper module to work with annotations (CVTerms)

Stores dictionaries which hold information the identifiers.org syntax, has functions to add CVTerms to different entities and parse CVTerms.
"""
from libsbml import CVTerm, BIOLOGICAL_QUALIFIER, BQB_IS

metabol_db_dict = {'BIGG': 'bigg.metabolite/',
                   'BRENDA': 'brenda/',
                   'CHEBI': 'chebi/',
                   'INCHI': 'inchi/',
                   'KEGG': 'kegg.compound/',
                   'METACYC': 'metacyc.compound/',
                   'MXNREF': 'metanetx.chemical/',
                   'SEED': 'seed.compound/',
                   'UPA': 'unipathway.compound/',
                   'HMDB': 'hmdb/',
                   'REACTOME': 'reactome/',
                   'BIOCYC': 'biocyc/',
                   'PUBCHEM': 'pubchem.compound/'}

reaction_db_dict = {'BIGG': 'bigg.reaction/',
                    'BRENDA': 'brenda/',
                    'KEGG': 'kegg.reaction/',
                    'METACYC': 'metacyc.reaction/',
                    'MXNREF': 'metanetx.reaction/',
                    'SEED': 'seed.reaction/',
                    'UPA': 'unipathway.reaction/',
                    'HMDB': 'hmdb/',
                    'REACTOME': 'reactome/',
                    'RHEA': 'rhea/',
                    'EC': 'ec-code/'}

gene_db_dict = {'NCBI': 'ncbiprotein/'}

pathway_db_dict = {'KEGG': 'kegg.pathway/'}


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
        all_ids.extend(ids)

    return all_ids
