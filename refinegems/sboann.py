#!/usr/bin/env python
"""Provides functions to automate the addition of SBO terms to the model

Script written by Elisabeth Fritze in her bachelor thesis.
Modified by Gwendolyn O. Döbel during her master thesis.
Commented by Famke Bäuerle and extended by Nantia Leonidou.

It is splitted into a lot of small functions which are all annotated, however when using it for SBO-Term annotation it only makes sense to run the "main" function: sbo_annotation(model_libsbml, database_user, database_name) if you want to continue with the model. The smaller functions might be useful if special information is needed for a reaction without the context of a bigger model or when the automated annotation fails for some reason.
"""

import re
import sqlite3
from libsbml import SpeciesReference, Compartment, Reaction, Species
from libsbml import Model as libModel
from refinegems.databases import PATH_TO_DB

__author__ = "Elisabeth Fritze, Gwendolyn O. Döbel, Famke Baeuerle and Nantia Leonidou"


def getCompartmentlessSpeciesId(speciesReference: SpeciesReference) -> str:
    """Determines wheter a species has compartment by its refernece

    Args:
        - speciesReference (SpeciesReference): Reference to species

    Returns:
        libsbml-species-id: id of species without compartment
    """
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    wasteStringLen = len(compartment) + 1
    return(speciesId[:-wasteStringLen])


def getCompartmentFromSpeciesRef(speciesReference: SpeciesReference) -> Compartment:
    """Extracts compartment from a species by its reference

    Args:
        - speciesReference (SpeciesReference): Reference to species

    Returns:
        Compartment: Compartment which the species lives in
    """
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    return compartment


def returnCompartment(id):
    """Helper to split compartment id"""
    return id[-1]


def getReactantIds(reac: Reaction) -> list[str]:
    """Extracts reactants (metabolites) of reaction

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list[str]: Reactants (metabolites) ids
    """
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(metabolite.getSpecies())
    return list


def getCompartmentlessReactantIds(reac: Reaction):
    """Extracts reactants which have no compartment information

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: reactants (metabolites) without compartments
    """
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list


def getProductIds(reac: Reaction):
    """Extracts products (metabolites) of reaction

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: products (metabolites) ids
    """
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(metabolite.getSpecies())
    return list


def getCompartmentlessProductIds(reac: Reaction):
    """Extracts products which have no compartment information

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: products (metabolites) without compartments
    """
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list


def getListOfMetabolites(reac: Reaction):
    """Extracts list of metabolites of the reaction

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: metabolites that are part of the reaction
    """
    list = []
    for reactant in reac.getListOfReactants():
        list.append(reactant)
    for product in reac.getListOfProducts():
        list.append(product)
    return list


def getMetaboliteIds(reac: Reaction):
    """Extracts list of metabolite ids of reaction

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: metabolite ids
    """
    return getReactantIds(reac) + getProductIds(reac)


def getCompartmentlessMetaboliteIds(reac: Reaction):
    """Extracts metabolites which have no compartment information

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: all metabolites which have no compartment
    """
    return getCompartmentlessReactantIds(
        reac) + getCompartmentlessProductIds(reac)


def getReactantCompartmentList(reac: Reaction):
    """Extracts compartments of reactants

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        set: compartment information of all reactants (metabolites)
    """
    compartments = []
    for metabolite in reac.getListOfReactants():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)


def getProductCompartmentList(reac: Reaction):
    """Extracts compartments of products

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        set: compartment information of all products (metabolites)
    """
    compartments = []
    for metabolite in reac.getListOfProducts():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)


def getCompartmentList(reac: Reaction):
    """Extracts compartments of metabolites

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        set: compartment information of all metabolites
    """
    metabolites = getListOfMetabolites(reac)
    compartments = []
    for metabolite in metabolites:
        compartments.append(getCompartmentFromSpeciesRef(metabolite))
    return set(compartments)


def getCompartmentDict(reac: Reaction):
    """sorts metabolites by compartment

    Args:
        - reac (Reaction): Reaction from sbml model

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


def moreThanTwoCompartmentTransport(reac: Reaction):
    """check if reaction traverses more than 2 compartments

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        bool: True if reaction traverses more than 2 compartments
    """
    return len(getCompartmentList(reac)) > 2


def isProtonTransport(reac: Reaction):
    """check if reaction is proton transport

    Args:
        - reac (Reaction): Reaction from sbml model

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


def soleProtonTransported(reac: Reaction):
    """check if reaction is transport powered by one H

    Args:
        - reac (Reaction): Reaction from sbml model

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


def getECNums(reac: Reaction):
    """Extracts EC-Code from the reaction annotations

    Args:
        - reac (Reaction): Reaction from sbml model

    Returns:
        list: all EC-Numbers of the reaction
    """
    lines = reac.getAnnotationString().split("\n")
    ECNums = []
    for line in lines:
        if 'ec-code' in line:
            try:
                ECNums.append(line.split('ec-code:')[1][:-3])
            except (IndexError):
                ECNums.append(line.split('ec-code/')[1][:-3])
    return ECNums


def splitTransportBiochem(reac: Reaction):
    """Tests if reaction traverses more than 1 compartment and set SBO Term

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getCompartmentList(reac)) > 1 and not soleProtonTransported(reac):
        reac.setSBOTerm('SBO:0000655')
    else:
        reac.setSBOTerm('SBO:0000176')


def checkSink(reac: Reaction):
    """Tests if reac is sink and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if "_SINK_" in reac.getId() or "_SK_" in reac.getId():
        reac.setSBOTerm('SBO:0000632')


def checkExchange(reac: Reaction):
    """Tests if reac is exchange and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if "_EX_" in reac.getId():
        reac.setSBOTerm('SBO:0000627')


def checkDemand(reac: Reaction):
    """Tests if reac is demand and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if "_DM_" in reac.getId():
        reac.setSBOTerm('SBO:0000628')


def checkBiomass(reac: Reaction):  # memote says growth is biomass
    """Tests if reac is biomass / growth and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    # Use regex to generalise check for growth/biomass reaction
    regex = 'growth|_*biomass\d*_*'
    if re.match(regex, reac.getId(), re.IGNORECASE):
        reac.setSBOTerm('SBO:0000629')


def checkPassiveTransport(reac: Reaction):
    """Tests if reac is passive transport and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    reactants = reac.getListOfReactants()
    products = reac.getListOfProducts()
    if len(reactants) == len(products) == 1:
        reactant = reactants[0].getSpecies()
        product = products[0].getSpecies()
        if returnCompartment(reactant) != returnCompartment(product):
            reac.setSBOTerm('SBO:0000658')


def checkActiveTransport(reac: Reaction):
    """Tests if reac is active transport (uses atp/pep) and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    reactantIds = []
    for metabolite in reac.getListOfReactants():
        reactantIds.append(metabolite.getSpecies())
    if 'M_atp_c' in reactantIds or 'M_pep_c' in reactantIds:
        reac.setSBOTerm('SBO:0000657')
        if reac.getReversible():
            print("Error, active reaction but reversible " + reac.getId())


def checkCoTransport(reac: Reaction):
    """Tests if reac is co-transport and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    reactants = reac.getListOfReactants()
    if len(reactants) > 1:
        reac.setSBOTerm('SBO:0000654')


def splitSymAntiPorter(reac: Reaction):
    """Tests if reac is sym- or antiporter and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getCompartmentList(reac)) > 2:
        pass
    elif 1 == len(getReactantCompartmentList(reac)) == len(getProductCompartmentList(reac)):
        reac.setSBOTerm('SBO:0000659')
    else:
        reac.setSBOTerm('SBO:0000660')


def checkPhosphorylation(reac: Reaction):
    """Tests if reac is phosphorylase / kinase and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    name = reac.getName()
    atpIsReactant = 'M_atp_c' in getReactantIds(reac)
    adpIsProduct = 'M_adp_c' in getProductIds(reac)
    if ' phosphorylase' in name or 'kinase' in name:
        reac.setSBOTerm('SBO:0000216')


def hasReactantPair(reac: Reaction, met1: Species, met2: Species) -> bool:
    """| Checks if a pair of metabolites is present in reaction
       | needed for special reactions like redox or deamination

    Args:
        - reac (Reaction): Reaction from sbml model
        - met1 (Species): metabolite 1 of metabolite pair
        - met2 (Species): metabolite 2 of metabolite pair

    Returns:
        bool: True if one of the metabolites is in reactants and the other in products
    """
    reactants = getCompartmentlessReactantIds(reac)
    products = getCompartmentlessProductIds(reac)
    return (met1 in reactants and met2 in products) or (
        met2 in reactants and met1 in products)


def checkRedox(reac: Reaction):
    """Tests if reac is redox and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
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


def checkGlycosylation(reac: Reaction):
    """Tests if reac is glycosylation and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    isGlycosylation = False
    isGlycosylation = isGlycosylation or hasReactantPair(
        reac, 'M_ppi', 'M_prpp')
    isGlycosylation = isGlycosylation or hasReactantPair(
        reac, 'M_udpglcur', 'M_udp')
    if isGlycosylation:
        reac.setSBOTerm('SBO:0000217')


def checkDecarbonylation(reac: Reaction):
    """Tests if reac is decarbonylation and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if not reac.getReversible():
        if 'M_co' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000400')


def checkDecarboxylation(reac: Reaction):
    """Tests if reac is decarboxylation and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if not reac.getReversible():
        if 'M_co2' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000399')


def checkDeamination(reac: Reaction):
    """Tests if reac is deamination and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if not reac.getReversible():
        waterAdded = 'M_h2o' in getCompartmentlessReactantIds(reac)
        nh4Removed = 'M_nh4' in getCompartmentlessProductIds(reac)
        if waterAdded and nh4Removed:
            reac.setSBOTerm('SBO:0000401')


def checkRedoxViaEC(reac: Reaction):
    """Tests if reac is redox by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('1'):
            reac.setSBOTerm('SBO:0000200')


def checkAcetylationViaEC(reac: Reaction):
    """Tests if reac is acetylation by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.3.1'):
            reac.setSBOTerm('SBO:0000215')


def checkGlycosylationViaEC(reac: Reaction):
    """Tests if reac is glycosylation by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.4'):
            reac.setSBOTerm('SBO:0000217')


def checkMethylationViaEC(reac: Reaction):
    """Tests if reac is methylation by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.1.1'):
            reac.setSBOTerm('SBO:0000214')


def checkTransaminationViaEC(reac: Reaction):
    """Tests if reac is transamination by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.6.1'):
            reac.setSBOTerm('SBO:0000403')


def checkDeaminationViaEC(reac: Reaction):
    """Tests if reac is deamination by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3.5.4'):
            reac.setSBOTerm('SBO:0000401')


def checkDecarboxylationViaEC(reac: Reaction):
    """Tests if reac is decarboxylation by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('4.1.1'):
            reac.setSBOTerm('SBO:0000399')


def checkIsomerisationViaEC(reac: Reaction):
    """Tests if reac is isomerisation by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('5'):
            reac.setSBOTerm('SBO:0000377')


def checkHydrolysisViaEC(reac: Reaction):
    """Tests if reac is hydrolysis by its EC-Code and sets SBO Term if true

    Args:
        - reac (Reaction): Reaction from sbml model
    """
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3'):
            reac.setSBOTerm('SBO:0000376')


def addSBOviaEC(reac: Reaction, cur):
    """Adds SBO terms based on EC numbers given in the annotations of a reactions

    Args:
        - reac (Reaction): Reaction from sbml model
        - cur (sqlite3.connect.cursor): Used to access the sqlite3 database
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


def addSBOfromDB(reac: Reaction, cur) -> bool:
    """Adds SBO term based on bigg id of a reaction

    Args:
        - reac (Reaction): Reaction from sbml model
        - cur (sqlite3.connect.cursor): Used to access the sqlite3 database

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

### functions below from Nantia ###

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


def addSBOforModel(model):
    model.setSBOTerm("SBO:0000624")


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

### end functions from Nantia ###

def sbo_annotation(model_libsbml: libModel) -> libModel:
    """Executes all steps to annotate SBO terms to a given model (former main function of original script by Elisabeth Fritze)

    Args:
        - model_libsbml (libModel): Model loaded with libsbml
      
    Returns:
        libModel: Modified model with SBO terms
    """
    open_con = sqlite3.connect(PATH_TO_DB)
    open_cur = open_con.cursor()

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
    
    ### functions from Nantia ###            
    # addSBOforMetabolites(model_libsbml)
    # addSBOforGenes(model_libsbml)
    # addSBOforModel(model_libsbml)
    # addSBOforGroups(model_libsbml)
    # addSBOforParameters(model_libsbml)
    # addSBOforCompartments(model_libsbml)
    
    open_cur.close()
    open_con.close()
    return model_libsbml