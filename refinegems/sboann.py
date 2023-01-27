#!/usr/bin/env python
"""Provides functions to automate the addition of SBO terms to the model

Script written by Elisabeth Fritze in her bachelor thesis.
Modified by Gwendolyn O. Gusak during her master thesis.
Commented by Famke Bäuerle and extended by Nantia Leonidou.

It is splitted into a lot of small functions which are all annotated, however when using it for SBO-Term annotation it only makes sense to run the "main" function: sbo_annotation_write(model_libsbml, database_user, database_name, new_filename) if you want to write the modified model to a SBML file or sbo_annotation(model_libsbml, database_user, database_name) if you want to continue with the model. The smaller functions might be useful if special information is needed for a reaction without the context of a bigger model or when the automated annotation fails for some reason.
"""

import re
import sqlite3
from sqlite3 import Error
from libsbml import *
from refinegems.load import write_to_file

__author__ = "Elisabeth Fritze"
__author__ = "Gwendolyn O. Gusak"


def is_valid_database(open_cur) -> bool:
   """
   Verifies if database has 2 tables with names 'bigg_to_sbo' & 'ec_to_sbo'
   
   Args:
        open_cur: sqlite3.connect.cursor() object of a database
 
   Returns:
        Boolean: True if all conditions (2 tables with names 'bigg_to_sbo' & 'ec_to_sbo') for database correct
   """

   # Fetches the table names as string tuples from the connected database
   open_cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
   tables = [string[0] for string in open_cur.fetchall()]

   return 'bigg_to_sbo' in tables and 'ec_to_sbo' in tables and len(tables) == 2


def initialise_SBO_database():
   """
   Initialises the SBO annotation database with 2 tables ('bigg_to_sbo' & 'ec_to_sbo')
   if file './data/sbo/sbo_database.db' is an incorrect database,
   otherwise correct database already exists
   
   Returns:
        sqlite3.connect() & sqlite3.connect.cursor() objects for SBO database
   """

   # Initialise empty connection
   con = None

   # Try to open connection & get cursor
   try:
      con = sqlite3.connect('./data/sbo/sbo_database.db')
      cursor = con.cursor()

      # If connected to database incorrect -> initialise correct one from file './data/sbo/sbo_database.sql'
      if not is_valid_database(cursor):
         with open('./data/sbo/sbo_database.sql') as schema:
            cursor.executescript(schema.read())
   except Error as e:
      print(e)
   finally:
      if con:
         return con, cursor


def getCompartmentlessSpeciesId(speciesReference):
    """determines wheter a species has compartment by its refernece

    Args:
        speciesReference (libsbml-speciesreference): reference to species

    Returns:
        libsbml-species-id: id of species without compartment
    """
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    wasteStringLen = len(compartment) + 1
    return(speciesId[:-wasteStringLen])


def getCompartmentFromSpeciesRef(speciesReference):
    """extracts compartment from a species by its reference

    Args:
        speciesReference (libsbml-speciesreference): reference to species

    Returns:
        libsbml-compartment: compartment which the species lives in
    """
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    return compartment


def returnCompartment(id):
    """helper to split compartment id"""
    return id[-1]


def getReactantIds(reac):
    """extracts reactants (metabolites) of reaction

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: reactants (metabolites) ids
    """
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(metabolite.getSpecies())
    return list


def getCompartmentlessReactantIds(reac):
    """extracts reactants which have no compartment information

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: reactants (metabolites) without compartments
    """
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list


def getProductIds(reac):
    """extracts products (metabolites) of reaction

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: products (metabolites) ids
    """
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(metabolite.getSpecies())
    return list


def getCompartmentlessProductIds(reac):
    """extracts products which have no compartment information

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: products (metabolites) without compartments
    """
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list


def getListOfMetabolites(reac):
    """extracts list of metabolites of the reaction

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: metabolites that are part of the reaction
    """
    list = []
    for reactant in reac.getListOfReactants():
        list.append(reactant)
    for product in reac.getListOfProducts():
        list.append(product)
    return list


def getMetaboliteIds(reac):
    """extracts list of metabolite ids of reaction

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: metabolite ids
    """
    return getReactantIds(reac) + getProductIds(reac)


def getCompartmentlessMetaboliteIds(reac):
    """extracts metabolites which have no compartment information

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: all metabolites which have no compartment
    """
    return getCompartmentlessReactantIds(
        reac) + getCompartmentlessProductIds(reac)


def getReactantCompartmentList(reac):
    """extracts compartments of reactants

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        set: compartment information of all reactants (metabolites)
    """
    compartments = []
    for metabolite in reac.getListOfReactants():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)


def getProductCompartmentList(reac):
    """extracts compartments of products

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        set: compartment information of all products (metabolites)
    """
    compartments = []
    for metabolite in reac.getListOfProducts():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)


def getCompartmentList(reac):
    """extracts compartments of metabolites

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        set: compartment information of all metabolites
    """
    metabolites = getListOfMetabolites(reac)
    compartments = []
    for metabolite in metabolites:
        compartments.append(getCompartmentFromSpeciesRef(metabolite))
    return set(compartments)


def getCompartmentDict(reac):
    """sorts metabolites by compartment

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        dict: compartment as key and metabolites as values
    """
    compartmentDict = {}
    for compartment in getCompartmentList(reac):
        compartmentDict[compartment] = []
        for metabolite in getListOfMetabolites(reac):
            if '_' + compartment in metabolite.getSpecies():
                compartmentDict[compartment].append(
                    getCompartmentlessSpeciesId(metabolite))
    return compartmentDict


def moreThanTwoCompartmentTransport(reac):
    """check if reaction traverses more than 2 compartments

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        bool: True if reaction traverses more than 2 compartments
    """
    return len(getCompartmentList(reac)) > 2


def isProtonTransport(reac):
    """check if reaction is proton transport

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        bool: True if reaction is proton transport
    """
    reactants = reac.getListOfReactants()
    products = reac.getListOfProducts()
    protonTransport = False
    if len(products) == len(reactants) == 1:
        if len(getCompartmentList(reac)) == 2:
            reactant = getCompartmentlessSpeciesId(reactants[0])
            product = getCompartmentlessSpeciesId(products[0])
            protonTransport = reactant == product == 'M_h'
    return protonTransport


def soleProtonTransported(reac):
    """check if reaction is transport powered by one H

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        bool: True if reaction is transport powered by one H
    """
    dict = getCompartmentDict(reac)
    soleProton = False
    for compartment in dict:
        if len(dict[compartment]) == 1:
            soleProton = soleProton or dict[compartment][0] == "M_h"
    return soleProton and not isProtonTransport(
        reac) and not moreThanTwoCompartmentTransport(reac)


def getECNums(reac):
    """extracts EC-Code from the reaction annotations

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model

    Returns:
        list: all EC-Numbers of the reaction
    """
    lines = reac.getAnnotationString().split("\n")
    ECNums = []
    for line in lines:
        if 'ec-code' in line:
            ECNums.append(line.split('ec-code/')[1][:-3])
    return ECNums


def splitTransportBiochem(reac):
    """test if reaction traverses more than 1 compartment and set SBO Term

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getCompartmentList(reac)) > 1 and not soleProtonTransported(reac):
        reac.setSBOTerm('SBO:0000655')
    else:
        reac.setSBOTerm('SBO:0000176')


def checkSink(reac):
    """tests if reac is sink and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if "_SINK_" in reac.getId() or "_SK_" in reac.getId():
        reac.setSBOTerm('SBO:0000632')


def checkExchange(reac):
    """tests if reac is exchange and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if "_EX_" in reac.getId():
        reac.setSBOTerm('SBO:0000627')


def checkDemand(reac):
    """tests if reac is demand and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if "_DM_" in reac.getId():
        reac.setSBOTerm('SBO:0000628')


def checkBiomass(reac):  # memote says growth is biomass
    """tests if reac is biomass / growth and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    # Use regex to generalise check for growth/biomass reaction
    regex = 'growth|_*biomass\d*_*'
    if re.match(regex, reac.getId(), re.IGNORECASE):
        reac.setSBOTerm('SBO:0000629')


def checkPassiveTransport(reac):
    """tests if reac is passive transport and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    reactants = reac.getListOfReactants()
    products = reac.getListOfProducts()
    if len(reactants) == len(products) == 1:
        reactant = reactants[0].getSpecies()
        product = products[0].getSpecies()
        if returnCompartment(reactant) != returnCompartment(product):
            reac.setSBOTerm('SBO:0000658')


def checkActiveTransport(reac):
    """tests if reac is active transport (uses atp/pep) and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    reactantIds = []
    for metabolite in reac.getListOfReactants():
        reactantIds.append(metabolite.getSpecies())
    if 'M_atp_c' in reactantIds or 'M_pep_c' in reactantIds:
        reac.setSBOTerm('SBO:0000657')
        if reac.getReversible():
            print("Error, active reaction but reversible " + reac.getId())


def checkCoTransport(reac):
    """tests if reac is co-transport and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    reactants = reac.getListOfReactants()
    if len(reactants) > 1:
        reac.setSBOTerm('SBO:0000654')


def splitSymAntiPorter(reac):
    """tests if reac is sym- or antiporter and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getCompartmentList(reac)) > 2:
        pass
    elif 1 == len(getReactantCompartmentList(reac)) == len(getProductCompartmentList(reac)):
        reac.setSBOTerm('SBO:0000659')
    else:
        reac.setSBOTerm('SBO:0000660')


def checkPhosphorylation(reac):
    """tests if reac is phosphorylase / kinase and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    name = reac.getName()
    atpIsReactant = 'M_atp_c' in getReactantIds(reac)
    adpIsProduct = 'M_adp_c' in getProductIds(reac)
    if ' phosphorylase' in name or 'kinase' in name:
        reac.setSBOTerm('SBO:0000216')


def hasReactantPair(reac, met1, met2):
    """checks if a pair of metabolites is present in reaction
       needed for special reactions like redox or deamination

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
        met1 (libsbml-metabolite): metabolite 1 of metabolite pair
        met2 (libsbml-metabolite): metabolite 2 of metabolite pair

    Returns:
        bool: True if one of the metabolites is in reactants and the other in products
    """
    reactants = getCompartmentlessReactantIds(reac)
    products = getCompartmentlessProductIds(reac)
    return (met1 in reactants and met2 in products) or (
        met2 in reactants and met1 in products)


def checkRedox(reac):
    """tests if reac is redox and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    isRedox = False
    isRedox = isRedox or hasReactantPair(reac, 'M_pyr', 'M_lac_L')
    isRedox = isRedox or hasReactantPair(reac, 'M_pyr', 'M_lac_D')
    isRedox = isRedox or hasReactantPair(reac, 'M_nad', 'M_nadh')
    isRedox = isRedox or hasReactantPair(reac, 'M_nadp', 'M_nadph')
    isRedox = isRedox or hasReactantPair(reac, 'M_fad', 'M_fadh2')
    isRedox = isRedox or hasReactantPair(reac, 'M_h20', 'M_h2o2')
    if isRedox:
        reac.setSBOTerm('SBO:0000200')


def checkGlycosylation(reac):
    """tests if reac is glycosylation and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    isGlycosylation = False
    isGlycosylation = isGlycosylation or hasReactantPair(
        reac, 'M_ppi', 'M_prpp')
    isGlycosylation = isGlycosylation or hasReactantPair(
        reac, 'M_udpglcur', 'M_udp')
    if isGlycosylation:
        reac.setSBOTerm('SBO:0000217')


def checkDecarbonylation(reac):
    """tests if reac is decarbonylation and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if not reac.getReversible():
        if 'M_co' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000400')


def checkDecarboxylation(reac):
    """tests if reac is decarboxylation and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if not reac.getReversible():
        if 'M_co2' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000399')


def checkDeamination(reac):
    """_tests if reac is deamination and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if not reac.getReversible():
        waterAdded = 'M_h2o' in getCompartmentlessReactantIds(reac)
        nh4Removed = 'M_nh4' in getCompartmentlessProductIds(reac)
        if waterAdded and nh4Removed:
            reac.setSBOTerm('SBO:0000401')


def checkRedoxViaEC(reac):
    """tests if reac is redox by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('1'):
            reac.setSBOTerm('SBO:0000200')


def checkAcetylationViaEC(reac):
    """tests if reac is acetylation by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.3.1'):
            reac.setSBOTerm('SBO:0000215')


def checkGlycosylationViaEC(reac):
    """tests if reac is glycosylation by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.4'):
            reac.setSBOTerm('SBO:0000217')


def checkMethylationViaEC(reac):
    """tests if reac is methylation by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.1.1'):
            reac.setSBOTerm('SBO:0000214')


def checkTransaminationViaEC(reac):
    """tests if reac is transamination by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.6.1'):
            reac.setSBOTerm('SBO:0000403')


def checkDeaminationViaEC(reac):
    """tests if reac is deamination by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3.5.4'):
            reac.setSBOTerm('SBO:0000401')


def checkDecarboxylationViaEC(reac):
    """tests if reac is decarboxylation by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('4.1.1'):
            reac.setSBOTerm('SBO:0000399')


def checkIsomerisationViaEC(reac):
    """tests if reac is isomerisation by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('5'):
            reac.setSBOTerm('SBO:0000377')


def checkHydrolysisViaEC(reac):
    """tests if reac is hydrolysis by its EC-Code and sets SBO Term if true

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3'):
            reac.setSBOTerm('SBO:0000376')


def addSBOviaEC(reac, cur):
    """Adds SBO terms based on EC numbers given in the annotations of a reactions

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
        cur (sqlite3.connect.cursor): used to access the sqlite3 database
    """
    if len(getECNums(reac)) == 1:
        ECnum = getECNums(reac)[0]
        splittedEC = ECnum.split(".")
        if len(splittedEC) == 4:
            ECpos1 = splittedEC[0]
            ECpos1to2 = ECpos1 + "." + splittedEC[1]
            ECpos1to3 = ECpos1to2 + "." + splittedEC[2]
            query4 = cur.execute("""SELECT sbo_term
                                     FROM ec_to_sbo 
                                    WHERE ecnum = ?""", [ECnum])
            result4 = cur.fetchone()
            if result4 is not None:
                sbo4 = result4[0]
                reac.setSBOTerm(sbo4)
            else:
                query3 = cur.execute("""SELECT sbo_term
                                        FROM ec_to_sbo 
                                        WHERE ecnum = ?""", [ECpos1to3])
                result3 = cur.fetchone()
                if result3 is not None:
                    sbo3 = result3[0]
                    reac.setSBOTerm(sbo3)
                else:
                    query2 = cur.execute("""SELECT sbo_term
                                            FROM ec_to_sbo t
                                            WHERE ecnum = ?""", [ECpos1to2])
                    result2 = cur.fetchone()
                    if result2 is not None:
                        sbo2 = result2[0]
                        reac.setSBOTerm(sbo2)
                    else:
                        query1 = cur.execute("""SELECT sbo_term
                                                FROM ec_to_sbo
                                                WHERE ecnum = ?""", [ECpos1])
                        result1 = cur.fetchone()
                        if result1 is not None:
                            sbo1 = result1[0]
                            reac.setSBOTerm(sbo1)


def addSBOfromDB(reac, cur):
    """Adds SBO term based on bigg id of a reaction

    Args:
        reac (libsbml-reaction): libsbml reaction from sbml model
        cur (sqlite3.connect.cursor): used to access the sqlite3 database

    Returns:
        bool: True if SBO Term was changed
    """
    reacid = reac.getId()
    query = cur.execute("""SELECT sbo_term 
                           FROM bigg_to_sbo 
                           WHERE bigg_reactionid = ?""", [reacid])
    result = cur.fetchone()
    if result is not None:
        sbo_term = result[0]
        reac.setSBOTerm(sbo_term)
        return True
    else:
        return False

### functions below from Nantia, not yet linked in main ###

def addSBOforMetabolites(model):
    # add metabolites SBO
    for met in model.species:
        met_id = met.getId()
        model.getSpecies(met_id).setSBOTerm("SBO:0000247")


def addSBOforGenes(model):
    # add genes SBO
    model_fbc = model.getPlugin("fbc")
    # if model has genes
    if model_fbc is not None:
        for gene in model_fbc.getListOfGeneProducts():
            gene.setSBOTerm("SBO:0000243")


def addSBOforModel(doc):
    doc.setSBOTerm("SBO:0000624")


def addSBOforGroups(model):
    mplugin = model.getPlugin("groups")
    # if groups are in model defined
    if mplugin is not None:
        for grp in mplugin.getListOfGroups():
            grp.setSBOTerm("SBO:0000633")


def addSBOforParameters(model):
    for param in model.getListOfParameters():
        # reaction bounds
        if 'R_' in param.getId():
            param.setSBOTerm("SBO:0000625")
        # default set bounds
        else:
            param.setSBOTerm("SBO:0000626")


def addSBOforCompartments(model):
    for cmp in model.getListOfCompartments():
        cmp.setSBOTerm("SBO:0000290")   

        
def handleMultipleECs(react, ECNums):
    # if no EC number annotated in model
    if len(ECNums) == 0:
        react.setSBOTerm('SBO:0000176')

    else:
        # store first digits of all annotated EC numbers
        lst = []
        for ec in ECNums:
            lst.append(ec.split(".")[0])

        # if ec numbers are from different enzyme classes, based on first digit
        # no ambiguous classification possible
        if len(set(lst)) > 1:
            react.setSBOTerm("SBO:0000176")  # metabolic rxn

        # if ec numbers are from the same enzyme classes,
        # assign parent SBO term based on first digit in EC number
        else:

            # Oxidoreductases
            if "1" in set(lst):
                react.setSBOTerm("SBO:0000200")
            # Transferase
            elif "2" in set(lst):
                react.setSBOTerm("SBO:0000402")
            # Hydrolases
            elif "3" in set(lst):
                react.setSBOTerm("SBO:0000376")
            # Lyases
            elif "4" in set(lst):
                react.setSBOTerm("SBO:0000211")
            # Isomerases
            elif "5" in set(lst):
                react.setSBOTerm("SBO:0000377")
            # Ligases, proper SBO is missing from graph --> use one for modification of covalent bonds
            elif "6" in set(lst):
                react.setSBOTerm("SBO:0000182")
            # Translocases
            elif "7" in set(lst):
                react.setSBOTerm("SBO:0000185")
            # Metabolic reactions
            else:
                react.setSBOTerm("SBO:0000176")

### end functions from Nantia ###

def sbo_annotation(model_libsbml):
    """executes all steps to annotate SBO terms to a given model
       (former main function of original script by Elisabeth Fritze)

    Args:
        model_libsbml (libsbml-model): model loaded with libsbml
      
    Returns:
        libsbml-model: modified model with SBO terms
    """
    open_con, open_cur = initialise_SBO_database()

    for reaction in model_libsbml.reactions:
        if not addSBOfromDB(reaction, open_cur):
            reaction.unsetSBOTerm()
            splitTransportBiochem(reaction)
            checkBiomass(reaction)
            checkSink(reaction)
            checkExchange(reaction)
            checkDemand(reaction)
            if reaction.getSBOTermID() == 'SBO:0000655': #transporter
                checkPassiveTransport(reaction)
                checkActiveTransport(reaction)
                if reaction.getSBOTermID() != 'SBO:0000657':
                    checkCoTransport(reaction)
                    if reaction.getSBOTermID() == 'SBO:0000654':
                        splitSymAntiPorter(reaction)
            if reaction.getSBOTermID() == 'SBO:0000176': #metabolic reaction
                addSBOviaEC(reaction, open_cur)
            if reaction.getSBOTermID() == 'SBO:0000176':
                checkRedox(reaction)
                checkGlycosylation(reaction)
                checkDecarbonylation(reaction)
                checkDecarboxylation(reaction)
                checkDeamination(reaction)
                checkPhosphorylation(reaction)
    
    open_cur.close()
    open_con.close()
    return model_libsbml
    

def sbo_annotation_write(model_libsbml, new_filename):
    """wrapper around sbo_annotation to include writing to a file

    Args:
        model_libsbml (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    new_model = sbo_annotation(model_libsbml)
    write_to_file(new_model, new_filename)
