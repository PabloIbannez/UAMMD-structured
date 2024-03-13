import os
import logging

import json
import jsbeautifier

import shutil

def add(args,
        COMPONENTS_PATH,DATA_PATH,
        UAMMD_STRUCTURED_PATH):

    logger = logging.getLogger("USCM")

    logger.info(f"Loading components file {COMPONENTS_PATH}/Components.json")
    try:
        with open(COMPONENTS_PATH+"/Components.json") as f:
            components = json.load(f)
    except:
        logger.error(f"Error loading {COMPONENTS_PATH}/Components.json")
        raise Exception(f"Error loading components file")

    logger.info(f"Loading data file {DATA_PATH}/Data.json")
    try:
        with open(DATA_PATH+"/Data.json") as f:
            data = json.load(f)
    except:
        logger.error(f"Error loading {DATA_PATH}/Data.json")
        raise Exception(f"Error loading data file")

    logger.info("\n")

    ADD = args.add

    #Check ADD has at least 4 elements
    if len(ADD) < 5:
        logger.error("Invalid number of arguments for add.")
        logger.error("Expected: component CLASS SUBCLASS TYPE SUBTYPE1 SUBTYPE2 ... or data CLASS NAME GETTER DTYPE DEFAULT(optional)")
        raise Exception("Invalid number of arguments for add.")

    ADD_TYPE = ADD[0]

    if ADD_TYPE == "component":

        logger.info("Adding a new component ...")

        if not args.file:
            logger.error("No file provided. Use --file to specify the file to add the component to.")
            raise Exception("No file provided")

        FILE = args.file[0]

        #Check if the file exists
        if not os.path.exists(FILE):
            logger.error(f"File {FILE} not found")
            raise Exception(f"File not found")

        CLASS    = ADD[1]
        SUBCLASS = ADD[2]
        TYPE     = ADD[3]
        SUBTYPES = ADD[4:]

        #Chekc if the combination CLASS SUBCLASS TYPE SUBTYPE is already in the components.json file
        for SUBTYPE in SUBTYPES:
            try:
                components[CLASS][SUBCLASS]
            except KeyError:
                continue

            for typ,subTyp,f in components[CLASS][SUBCLASS]:
                if typ == TYPE and subTyp == SUBTYPE:
                    logger.error(f"Component {CLASS},{SUBCLASS},{TYPE},{SUBTYPE} already exists")
                    raise Exception(f"Component already exists")

        for SUBTYPE in SUBTYPES:
            logger.info(f"Adding component {CLASS} {SUBCLASS} {TYPE} {SUBTYPE} from file {FILE}")

            #Check if CLASS and SUBCLASS are in the components dictionary
            if CLASS not in components:
                logger.error(f"Type {CLASS} not found in components dictionary")
                raise Exception(f"Type not found in components dictionary")

            if SUBCLASS not in components[CLASS]:
                #Add it to the dictionary
                components[CLASS][SUBCLASS] = []

            #Add the component to the dictionary
            components[CLASS][SUBCLASS].append([TYPE,SUBTYPE,FILE])

            #Update the components.json file. Using jsbeautifier to make it more readable
            with open(COMPONENTS_PATH+"/Components.json", "w") as f:
                f.write(jsbeautifier.beautify(json.dumps(components)))

            ########################################################

            fileFolder = UAMMD_STRUCTURED_PATH+"/"+CLASS+"/"+SUBCLASS+"/"+TYPE
            filePath   = fileFolder+"/"+FILE

            #Check if the folder exists. If not, create it
            if not os.path.exists(fileFolder):
                os.makedirs(fileFolder)

            #Move the file to the folder
            shutil.move(FILE,filePath)

    elif ADD_TYPE == "data":

        logger.info("Adding a new data ...")

        CLASS   = ADD[1]
        NAME    = ADD[2]
        GETTER  = ADD[3]
        DTYPE   = ADD[4]
        DEFAULT = None

        if len(ADD) == 6:
            DEFAULT = ADD[5]

        #Check if CLASS is in the data dictionary
        if CLASS not in data:
            #Error
            logger.error(f"Type {CLASS} not found in data dictionary")
            raise Exception(f"Type not found in data dictionary")

        #Check if NAME is already in the data dictionary
        for dataInfo in data[CLASS]:
            if dataInfo[0] == NAME:
                logger.error(f"Data {CLASS} {NAME} already exists")
                raise Exception(f"Data already exists")

        logger.info(f"Adding data {CLASS} {NAME} {GETTER} {DTYPE} {DEFAULT}")

        if DEFAULT is None:
            data[CLASS].append([NAME,GETTER,DTYPE])
        else:
            data[CLASS].append([NAME,GETTER,DTYPE,DEFAULT])

        with open(DATA_PATH+"/Data.json", "w") as f:
            f.write(jsbeautifier.beautify(json.dumps(data)))

    else:
        logger.error(f"Invalid add type {ADD_TYPE}")
        raise Exception(f"Invalid add type")
