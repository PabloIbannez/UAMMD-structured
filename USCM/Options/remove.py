import os
import logging

import json
import jsbeautifier

def remove(args,
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

    REMOVE = args.remove

    #Check if REMOVE has 3 elements
    if len(REMOVE) < 3:
        logger.error("Invalid number of arguments for remove.")
        logger.error("Expected: component CLASS SUBCLASS TYPE SUBTYPE ... or data CLASS NAME")
        raise Exception("Invalid number of arguments for remove.")

    REMOVE_TYPE = REMOVE[0]

    if REMOVE_TYPE == "component":

        logger.info("Removing a component ...\n")

        if len(REMOVE) < 5:
            logger.error("Invalid number of arguments for remove component.")
            logger.error("Expected: component CLASS SUBCLASS TYPE SUBTYPE ...")
            raise Exception("Invalid number of arguments for remove component.")

        CLASS    = REMOVE[1]
        SUBCLASS = REMOVE[2]
        TYPE   = REMOVE[3]
        SUBTYPE= REMOVE[4]

        #Check if the combination CLASS SUBCLASS TYPE SUBTYPE is in the components.json file
        try:
            components[CLASS][SUBCLASS]
        except KeyError:
            logger.error(f"Component {CLASS},{SUBCLASS},{TYPE},{SUBTYPE} not found")
            raise Exception(f"Component not found")

        found = False
        for typ,subTyp,f in components[CLASS][SUBCLASS]:
            if typ == TYPE and subTyp == SUBTYPE:
                found = True
                break

        if not found:
            logger.error(f"Component {CLASS},{SUBCLASS},{TYPE},{SUBTYPE} not found")
            raise Exception(f"Component not found")

        ####

        for typ,subTyp,f in components[CLASS][SUBCLASS]:
            if typ == TYPE and subTyp == SUBTYPE:
                logger.info(f"Removing component {CLASS} {SUBCLASS} {TYPE} {SUBTYPE}")
                components[CLASS][SUBCLASS].remove([typ,subTyp,f])
                break

        fileName = f
        fileAppearTwice = False

        for typ,subTyp,f in components[CLASS][SUBCLASS]:
            if f == fileName:
                fileAppearTwice = True
                break

        if not fileAppearTwice:
            #Remove the file
            fileFolder = UAMMD_STRUCTURED_PATH+"/"+CLASS+"/"+SUBCLASS+"/"+TYPE
            filePath   = fileFolder+"/"+fileName

            if os.path.exists(filePath):
                logger.info(f"Removing file {filePath}")
                os.remove(filePath)
            else:
                logger.warning(f"File {filePath} not found")
        else:
            logger.info(f"File {fileName} is used by other components (not removed)")

        if len(components[CLASS][SUBCLASS]) == 0:
            #Remove SUBCLASS entry
            logger.info(f"Removing subtype {SUBCLASS}")
            del components[CLASS][SUBCLASS]

        #Update the components.json file. Using jsbeautifier to make it more readable
        with open(COMPONENTS_PATH+"/Components.json", "w") as f:
            f.write(jsbeautifier.beautify(json.dumps(components)))

    elif REMOVE_TYPE == "data":

        logger.info("Removing data ...\n")

        CLASS = REMOVE[1]
        NAME = REMOVE[2]

        #Check if CLASS is in the data dictionary
        if CLASS not in data:
            #Error
            logger.error(f"Type {CLASS} not found in data dictionary")
            raise Exception(f"Type not found in data dictionary")

        #Check if NAME is already in the data dictionary
        found = False
        for dataInfo in data[CLASS]:
            if dataInfo[0] == NAME:
                found = True
                break

        if not found:
            logger.error(f"Data {CLASS} {NAME} not found")
            raise Exception(f"Data not found")

        for dataInfo in data[CLASS]:
            if dataInfo[0] == NAME:
                logger.info(f"Removing data {CLASS} {NAME}")
                data[CLASS].remove(dataInfo)
                break

        #Update the data.json file. Using jsbeautifier to make it more readable
        with open(DATA_PATH+"/Data.json", "w") as f:
            f.write(jsbeautifier.beautify(json.dumps(data)))

    else:
        logger.error(f"Invalid remove type {REMOVE_TYPE}")
        raise Exception(f"Invalid remove type")

