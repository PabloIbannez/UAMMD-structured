import os
import sys

import logging

import json

from colorama import init as colorama_init
from colorama import Fore
from colorama import Style

#Import generators
from .DataGeneration.particleDataGenerator import generateParticleData
from .DataGeneration.stateLoaderGenerator  import generateStateLoader

#Global data
from .DataGeneration.GlobalData.handlerGenerator import generateHandler
from .ComponentsGeneration.GlobalData.genericLoaderGenerator  import generateGenericLoader

from .ComponentsGeneration.Integrator.genericIntegratorLoaderGenerator import generateGenericIntegratorLoader

from .ComponentsGeneration.DataStructures.VerletConditionalListSet.verletConditionalListSetLoadersGenerator import generateVerletConditionalListSetLoaders
from .ComponentsGeneration.SimulationStep.genericSimulationStepLoadersGenerator                             import generateGenericSimulationStepLoader

from .ComponentsGeneration.Interactor.Patches.genericPatchesPotentialLoaderGenerator import generateGenericPatchesPotentialLoader
from .ComponentsGeneration.Interactor.genericPotentialLoaderGenerator import generateGenericPotentialLoader

###################################################################

def checkFoldersAndFiles(componentsPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    with open(componentsPath) as f:
        components = json.load(f)

    #Check if component folders exists
    folders2Check = ()
    filesInFolder = {}
    for CLASS in components:
        for SUBCLASS in components[CLASS]:
            classSet = []
            for TYPE,_,_ in components[CLASS][SUBCLASS]:
                classSet.append(TYPE)
            classSet = set(classSet)

            files    = {}
            for TYPE in classSet:
                files[TYPE] = []
            for TYPE,_,FILE in components[CLASS][SUBCLASS]:
                files[TYPE].append(FILE)

            for TYPE in classSet:
                folder = os.path.join(uammdStructuredBaseFolder,CLASS,SUBCLASS,TYPE)
                filesInFolder[folder] = files[TYPE]
                folders2Check += (folder,)

    allFoldersExist = True
    missedFolders = []
    for folder in folders2Check:
        if not os.path.exists(folder):
            if folder not in missedFolders:
                missedFolders.append(folder)
                logger.error(f"Folder \"{folder}\" does not exist")
            allFoldersExist = False

    if allFoldersExist:
        allFilesExist = True
        missedFiles = []
        for folder in folders2Check:
            #Check if files in folder
            for file in filesInFolder[folder]:
                if not os.path.exists(folder + "/" + file):
                    if folder + "/" + file not in missedFiles:
                        missedFiles.append(folder + "/" + file)
                        logger.error(f"File \"{folder + '/' + file}\" does not exist")
                    allFilesExist = False
    else:
        allFilesExist = False

    return allFoldersExist, allFilesExist


def checkIncludeInFile(include,file):

    logger = logging.getLogger("USCM")

    #Check if include is a list, if not, make it a list
    if not isinstance(include,list):
        include = [include]

    #Check filePath line is present in uammddStructuredBaseIncludeFolder
    for inc in include:
        with open(file,"r") as fuammd:
            isPresent = False
            for line in fuammd:
                if inc in line:
                    isPresent = True
                    break
            if not isPresent:
                logger.warning(f"\"{file}\" does not contain the line \"{inc}\"")


def generateIncludeFile(componentsPath,CLASS,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    with open(componentsPath,"r") as f:
        components = json.load(f)

    includeFiles = []
    #Iterate over all components in components
    for SUBCLASS in sorted(components[CLASS]):

        TYPES = []
        for TYPE,_,_ in sorted(components[CLASS][SUBCLASS]):
            TYPES.append(TYPE)
        TYPES = list(set(TYPES))

        filePath = f"{CLASS}/{SUBCLASS}/{SUBCLASS}Includers.cuh"

        logger.debug(f"Generating {filePath}")
        with open(uammdStructuredBaseFolder+"/"+filePath,"w") as fincl:
            fincl.write(f"#ifndef __{SUBCLASS.upper()}_INCLUDERS_CUH__\n")
            fincl.write(f"#define __{SUBCLASS.upper()}_INCLUDERS_CUH__\n\n")

            includedFiles = []

            #Check if the file SUBCLASS.cuh exists in CLASS/SUBCLASS
            if os.path.isfile(uammdStructuredBaseFolder+"/"+CLASS+"/"+SUBCLASS+"/"+SUBCLASS+".cuh"):
                logger.debug(Fore.YELLOW + f"Adding " + Fore.RESET + f"common file (SUBCLASS) {SUBCLASS}.cuh to {filePath}")
                file = f"\"{SUBCLASS}.cuh\""
                if file not in includedFiles:
                    includedFiles.append(file)
                    fincl.write(f"#include {file}\n")

            for TYPE in sorted(TYPES):

                #Check if the file TYPE.cuh exists in CLASS/SUBCLASS/TYPE/
                if os.path.isfile(uammdStructuredBaseFolder+"/"+CLASS+"/"+SUBCLASS+"/"+TYPE+"/"+TYPE+".cuh"):
                    logger.debug(Fore.YELLOW + f"Adding " + Fore.RESET + f"common file (TYPE) {TYPE}.cuh to {filePath}")
                    file = f"\"{TYPE}/{TYPE}.cuh\""
                    if file not in includedFiles:
                        includedFiles.append(file)
                        fincl.write(f"#include {file}\n")

                for COMP_TYPE,COMP_SUBTYPE,COMP_FILE in sorted(components[CLASS][SUBCLASS]):
                    if COMP_TYPE == TYPE:
                        file = f"\"{COMP_TYPE}/{COMP_FILE}\""
                        if file not in includedFiles:
                            includedFiles.append(file)
                            fincl.write(f"#include {file}\n")

            fincl.write(f"\n#endif //__{SUBCLASS.upper()}_INCLUDERS_CUH__\n")

            #Add #include to filePath
            filePath = "#include\""+filePath+"\""
            includeFiles.append(filePath)

    return includeFiles

###################################################################

def generate(args,
             COMPONENTS_PATH,DATA_PATH,
             UAMMD_STRUCTURED_PATH,UAMMD_STRUCTURED_INCLUDE):

    logger = logging.getLogger("USCM")

    logger.info(f"Loading components file {COMPONENTS_PATH}/Components.json")
    try:
        with open(COMPONENTS_PATH+"/Components.json") as f:
            components = json.load(f)
    except:
        logger.error(f"Error loading {COMPONENTS_PATH}/Components.json")
        raise Exception(f"Error loading components file")

    logger.info(f"Loading data file {DATA_PATH}/Data.json\n")
    try:
        with open(DATA_PATH+"/Data.json") as f:
            data = json.load(f)
    except:
        logger.error(f"Error loading {DATA_PATH}/Data.json")
        raise Exception(f"Error loading data file")

    #Check if folders and files exist
    allFoldersExist, allFilesExist = checkFoldersAndFiles(COMPONENTS_PATH+"/Components.json",UAMMD_STRUCTURED_PATH)
    if not allFoldersExist:
        logger.error(Fore.RED + "Some folders are missing. Please check the output of the script."+ Fore.RESET)
        raise Exception("Some folders are missing. Please check the output of the script.")
    if not allFilesExist:
        logger.error(Fore.RED + "Some files are missing. Please check the output of the script."+ Fore.RESET)
        raise Exception("Some files are missing. Please check the output of the script.")

    if allFoldersExist and allFilesExist:
        logger.debug(Fore.GREEN + "All folders and files exist."+ Fore.RESET+"\n")

    logger.info("Generating data ...")

    ##Generate global data handlers

    ##Generate units
    generateHandler("units",DATA_PATH+"/Data.json","Units",UAMMD_STRUCTURED_PATH+"/GlobalData/Units/","UnitsHandler.cuh",
                    additionalIncludes=[],
                    genGetter=True,genSetter=False)

    ##Generate fundamental
    generateHandler("fundamental",DATA_PATH+"/Data.json","Fundamental",UAMMD_STRUCTURED_PATH+"/GlobalData/Fundamental/","FundamentalHandler.cuh",
                    additionalIncludes=[],
                    genGetter=True,genSetter=True)

    ##Generate ensemble
    generateHandler("ensemble",DATA_PATH+"/Data.json","Ensemble",UAMMD_STRUCTURED_PATH+"/GlobalData/Ensemble/","EnsembleHandler.cuh",
                    additionalIncludes=["utils/Box.cuh"],
                    genGetter=True,genSetter=True)



    ##Generate data

    UAMMD_PARTICLE_DATA = [
        ["id",     "Id",     "int"],
        ["mass",   "Mass",   "real"],
        ["force",  "Force",  "real4"],
        ["virial", "Virial", "real"],
        ["energy", "Energy", "real"],
        ["vel",    "Vel",    "real3"],
        ["radius", "Radius", "real"],
        ["charge", "Charge", "real"],
        ["torque", "Torque", "real4"],
        ["angVel", "AngVel", "real4"],
        ["dir",    "Dir",    "real4"]
    ]
    generateParticleData(DATA_PATH+"/Data.json","ParticleData",UAMMD_STRUCTURED_PATH+"/../","preamble.h")

    #Check all data in state is also present in particle data
    for state in data["State"]:
        found = False
        for particleData in data["ParticleData"]+UAMMD_PARTICLE_DATA:
            if state[1] == particleData[1] and state[2] == particleData[2]:
                found = True
                break
        if not found:
            logger.error(f"State {state} not found in ParticleData. Please add it before generating the code.")
            raise Exception(f"State not found in ParticleData")

    #Check no default data is present in particle data and state
    for d in data["ParticleData"]+data["State"]:
        if len(d) > 3:
            logger.error(f"Default data {d} found in ParticleData/State.")
            raise Exception(f"Default data found in ParticleData/State")

    generateStateLoader(DATA_PATH+"/Data.json","State",UAMMD_STRUCTURED_PATH+"/ParticleData/","StateLoader.cuh")

    logger.info(Fore.GREEN + "Generation finished successfully."+ Fore.RESET+"\n")

    logger.info("Generating components ...")

    #Global data
    checkIncludeInFile(generateIncludeFile(COMPONENTS_PATH+"/Components.json",
                                           "GlobalData",
                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericLoader("types",COMPONENTS_PATH+"/Components.json",
                                             ["GlobalData","Types"],
                                             UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericLoader("units",COMPONENTS_PATH+"/Components.json",
                                             ["GlobalData","Units"],
                                             UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericLoader("fundamental",COMPONENTS_PATH+"/Components.json",
                                             ["GlobalData","Fundamental"],
                                             UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericLoader("ensemble",COMPONENTS_PATH+"/Components.json",
                                             ["GlobalData","Ensemble"],
                                             UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    #Integrator

    checkIncludeInFile(generateIncludeFile(COMPONENTS_PATH+"/Components.json",
                                           "Integrator",
                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericIntegratorLoader(COMPONENTS_PATH+"/Components.json",
                                                       ["Integrator"],
                                                       UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    #DataStructures

    checkIncludeInFile(generateIncludeFile(COMPONENTS_PATH+"/Components.json",
                                           "DataStructures",
                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    #checkIncludeInFile(generateVerletConditionalListSetLoaders(COMPONENTS_PATH+"/Components.json",
    #                                                           ["DataStructures","VerletConditionalListSet"],
    #                                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    #SimulationStep
    checkIncludeInFile(generateIncludeFile(COMPONENTS_PATH+"/Components.json",
                                           "SimulationStep",
                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)
    checkIncludeInFile(generateGenericSimulationStepLoader(COMPONENTS_PATH+"/Components.json",
                                                           ["SimulationStep"],
                                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    #Interactors
    checkIncludeInFile(generateIncludeFile(COMPONENTS_PATH+"/Components.json",
                                           "Interactor",
                                           UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)


    checkIncludeInFile(generateGenericPatchesPotentialLoader(COMPONENTS_PATH+"/Components.json",
                                                             ["Interactor","Patches"],
                                                             UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    checkIncludeInFile(generateGenericPotentialLoader(COMPONENTS_PATH+"/Components.json",
                                                      ["Interactor"],
                                                      UAMMD_STRUCTURED_PATH),UAMMD_STRUCTURED_INCLUDE)

    logger.info(Fore.GREEN + "Generation finished successfully."+ Fore.RESET+"\n")

