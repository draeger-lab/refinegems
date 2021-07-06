from libsbml import *
import psycopg2

__author__ = "Elisabeth Fritze"

def getCompartmentlessSpeciesId(speciesReference):
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    wasteStringLen = len(compartment) + 1
    return(speciesId[:-wasteStringLen])

def getCompartmentFromSpeciesRef(speciesReference):
    speciesId = speciesReference.getSpecies()
    species = speciesReference.getModel().getSpecies(speciesId)
    compartment = species.getCompartment()
    return compartment

def returnCompartment(id):
    return id[-1]

def getReactantIds(reac):
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(metabolite.getSpecies())
    return list


def getCompartmentlessReactantIds(reac):
    list = []
    for metabolite in reac.getListOfReactants():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list

def getProductIds(reac):
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(metabolite.getSpecies())
    return list

def getCompartmentlessProductIds(reac):
    list = []
    for metabolite in reac.getListOfProducts():
        list.append(getCompartmentlessSpeciesId(metabolite))
    return list

def getListOfMetabolites(reac):
    list = []
    for reactant in reac.getListOfReactants():
        list.append(reactant)
    for product in reac.getListOfProducts():
        list.append(product)
    return list

def getMetaboliteIds(reac):
    return getReactantIds(reac) + getProductIds(reac)

def getCompartmentlessMetaboliteIds(reac):
    return getCompartmentlessReactantIds(reac) + getCompartmentlessProductIds(reac)

def getReactantCompartmentList(reac):
    compartments = []
    for metabolite in reac.getListOfReactants():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)

def getProductCompartmentList(reac):
    compartments = []
    for metabolite in reac.getListOfProducts():
        compartment = getCompartmentFromSpeciesRef(metabolite)
        compartments.append(compartment)
    return set(compartments)

def getCompartmentList(reac):
    metabolites = getListOfMetabolites(reac)
    compartments = []
    for metabolite in metabolites:
        compartments.append(getCompartmentFromSpeciesRef(metabolite))
    return set(compartments)


def getCompartmentDict(reac):
    compartmentDict = {}
    for compartment in getCompartmentList(reac):
        compartmentDict[compartment]=[]
        for metabolite in getListOfMetabolites(reac):
            if '_' + compartment in metabolite.getSpecies():
                compartmentDict[compartment].append(getCompartmentlessSpeciesId(metabolite))
    return compartmentDict

def moreThanTwoCompartmentTransport(reac):
    return len(getCompartmentList(reac)) > 2

def isProtonTransport(reac):
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
    dict = getCompartmentDict(reac)
    soleProton = False
    for compartment in dict:
        if len(dict[compartment]) == 1:
            soleProton = soleProton or dict[compartment][0] == "M_h"
    return soleProton and not isProtonTransport(reac) and not moreThanTwoCompartmentTransport(reac)

def getECNums(reac):
    lines = reac.getAnnotationString().split("\n")
    ECNums = []
    for line in lines:
        if 'ec-code' in line:
            ECNums.append(line.split('ec-code/')[1][:-3])
    return ECNums

def splitTransportBiochem(reac):
    if len(getCompartmentList(reac)) > 1 and not soleProtonTransported(reac):
        reac.setSBOTerm('SBO:0000655')
    else:
        reac.setSBOTerm('SBO:0000176')


def checkSink(reac):
    if "_SINK_" in reac.getId() or "_SK_" in reac.getId():
        reac.setSBOTerm('SBO:0000632')

def checkExchange(reac):
    if "_EX_" in reac.getId():
        reac.setSBOTerm('SBO:0000627')
def checkDemand(reac):
    if "_DM_" in reac.getId():
        reac.setSBOTerm('SBO:0000628')
def checkBiomass(reac):
    if '_BIOMASS_' in reac.getId():
        reac.setSBOTerm('SBO:0000629')

def checkPassiveTransport(reac):
    reactants = reac.getListOfReactants()
    products = reac.getListOfProducts()
    if len(reactants) == len(products) == 1:
        reactant = reactants[0].getSpecies()
        product = products[0].getSpecies()
        if returnCompartment(reactant) != returnCompartment(product):
            reac.setSBOTerm('SBO:0000658')

def checkActiveTransport(reac):
    reactantIds = []
    for metabolite in reac.getListOfReactants():
        reactantIds.append(metabolite.getSpecies())
    if 'M_atp_c' in reactantIds or 'M_pep_c' in reactantIds:
        reac.setSBOTerm('SBO:0000657')
        if reac.getReversible():
            print("Error, active reaction but reversible " + reac.getId())

def checkCoTransport(reac):
    reactants = reac.getListOfReactants()
    if len(reactants) > 1:
        reac.setSBOTerm('SBO:0000654')

def splitSymAntiPorter(reac):
    if len(getCompartmentList(reac)) > 2:
        pass
    elif 1 == len(getReactantCompartmentList(reac)) == len(getProductCompartmentList(reac)):
        reac.setSBOTerm('SBO:0000659')
    else:
        reac.setSBOTerm('SBO:0000660')


def checkPhosphorylation(reac):
    name = reac.getName()
    atpIsReactant = 'M_atp_c' in getReactantIds(reac)
    adpIsProduct = 'M_adp_c' in getProductIds(reac)
    if ' phosphorylase' in name or 'kinase' in name:
        reac.setSBOTerm('SBO:0000216')

def hasReactantPair(reaction, met1, met2):
    reactants = getCompartmentlessReactantIds(reaction)
    products = getCompartmentlessProductIds(reaction)
    return (met1 in reactants and met2 in products) or (met2 in reactants and met1 in products)

def checkRedox(reac):
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
    isGlycosylation = False
    isGlycosylation = isGlycosylation or hasReactantPair(reac, 'M_ppi', 'M_prpp')
    isGlycosylation = isGlycosylation or hasReactantPair(reac, 'M_udpglcur', 'M_udp')
    if isGlycosylation:
        reac.setSBOTerm('SBO:0000217')

def checkDecarbonylation(reac):
    if not reac.getReversible():
        if 'M_co' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000400')

def checkDecarboxylation(reac):
    if not reac.getReversible():
        if 'M_co2' in getCompartmentlessProductIds(reac):
            reac.setSBOTerm('SBO:0000399')

def checkDeamination(reac):
    if not reac.getReversible():
        waterAdded = 'M_h2o' in getCompartmentlessReactantIds(reac)
        nh4Removed = 'M_nh4' in getCompartmentlessProductIds(reac)
        if waterAdded and nh4Removed:
            reac.setSBOTerm('SBO:0000401')

def checkRedoxViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('1'):
            reac.setSBOTerm('SBO:0000200')

def checkAcetylationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.3.1'):
            reac.setSBOTerm('SBO:0000215')

def checkGlycosylationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.4'):
            reac.setSBOTerm('SBO:0000217')

def checkMethylationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.1.1'):
            reac.setSBOTerm('SBO:0000214')


def checkTransaminationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('2.6.1'):
            reac.setSBOTerm('SBO:0000403')

def checkDeaminationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3.5.4'):
            reac.setSBOTerm('SBO:0000401')

def checkDecarboxylationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('4.1.1'):
            reac.setSBOTerm('SBO:0000399')

def checkIsomerisationViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('5'):
            reac.setSBOTerm('SBO:0000377')

def checkHydrolysisViaEC(reac):
    if len(getECNums(reac)) == 1:
        if getECNums(reac)[0].startswith('3'):
            reac.setSBOTerm('SBO:0000376')

def addSBOviaEC(reac, cur):
    if len(getECNums(reac)) == 1:
        ECnum = getECNums(reac)[0]
        splittedEC = ECnum.split(".")
        if len(splittedEC) == 4:
            ECpos1 = splittedEC[0]
            ECpos1to2 = ECpos1 + "." + splittedEC[1]
            ECpos1to3 = ECpos1to2 + "." + splittedEC[2]
            query4 = cur.execute("""SELECT t.sbo_term
                                     FROM ec_to_sbo t
                                    WHERE t.ecnum = %s""", (ECnum,))
            result4 = cur.fetchone()
            if result4 != None:
                sbo4 = result4[0]
                reac.setSBOTerm(sbo4)
            else:
                query3 = cur.execute("""SELECT t.sbo_term
                                          FROM ec_to_sbo t
                                         WHERE t.ecnum = %s""", (ECpos1to3,))
                result3 = cur.fetchone()
                if result3 != None:
                    sbo3 = result3[0]
                    reac.setSBOTerm(sbo3)
                else:
                    query2 = cur.execute("""SELECT t.sbo_term
                                              FROM ec_to_sbo t
                                             WHERE t.ecnum = %s""", (ECpos1to2,))
                    result2 = cur.fetchone()
                    if result2 != None:
                        sbo2 = result2[0]
                        reac.setSBOTerm(sbo2)
                    else:
                        query1 = cur.execute("""SELECT t.sbo_term
                                                  FROM ec_to_sbo t
                                                 WHERE t.ecnum = %s""", (ECpos1,))
                        result1 = cur.fetchone()
                        if result1 != None:
                            sbo1 = result1[0]
                            reac.setSBOTerm(sbo1)
def addSBOfromDB(reac, cur):
    reacid = reac.getId()  
    query = cur.execute("SELECT t.sbo_term \
                       FROM bigg_to_sbo t \
                      WHERE  t.bigg_reactionid = %s", (reacid,))
    result = cur.fetchone()
    if result != None:
        sbo_term = result[0]
        reac.setSBOTerm(sbo_term)
        return True
    else:
        return False
    
def write_to_file(model, new_filename):
    """Writes modified model to new file

    Args:
        model (libsbml-model): model loaded with libsbml
        new_filename (Str): filename for modified model
    """
    new_document = model.getSBMLDocument()
    writeSBMLToFile(new_document, new_filename)
    print("Model with SBO Annotations written to " + new_filename)

def sbo_annotation(model_libsbml, database_user, database_name, new_filename):
    """executes all steps to annotate sbo terms to a given model
       (former main function of original script by Elisabeth Fritze)
    
    Args:
        model_libsbml (libsbml-model): model loaded with libsbml
        database_user (Str): username for postgresql database
        database_name (Str): name of database in which create_dbs.sql was imported
        new_filename (Str): filename for modified model
    """
    conn = psycopg2.connect(dbname=database_name, user=database_user)
    cur = conn.cursor()

    conn.autocommit = True
    for reaction in model_libsbml.reactions:
        if not addSBOfromDB(reaction, cur):
            reaction.unsetSBOTerm()
            splitTransportBiochem(reaction)
            checkBiomass(reaction)
            checkSink(reaction)
            checkExchange(reaction)
            checkDemand(reaction)
            if reaction.getSBOTermID() == 'SBO:0000655':
                checkPassiveTransport(reaction)
                checkActiveTransport(reaction)
                if reaction.getSBOTermID() != 'SBO:0000657':
                    checkCoTransport(reaction)
                    if reaction.getSBOTermID() == 'SBO:0000654':
                        splitSymAntiPorter(reaction)
            if reaction.getSBOTermID() == 'SBO:0000176':
                addSBOviaEC(reaction, cur)
            if reaction.getSBOTermID() == 'SBO:0000176':
                checkRedox(reaction)
                checkGlycosylation(reaction)
                checkDecarbonylation(reaction)
                checkDecarboxylation(reaction)
                checkDeamination(reaction)
                checkPhosphorylation(reaction)
    
    write_to_file(model_libsbml, new_filename)
    cur.close()
    conn.close()
