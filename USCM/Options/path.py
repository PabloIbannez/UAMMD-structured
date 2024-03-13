import os
import logging

import json

def path(args,
         COMPONENTS_PATH,
         UAMMD_STRUCTURED_PATH):

    logger = logging.getLogger("USCM")

    logger.info("Getting the path of a component ...")

    logger.info(f"Loading components file {COMPONENTS_PATH}/Components.json")
    try:
        with open(COMPONENTS_PATH+"/Components.json") as f:
            components = json.load(f)
    except:
        logger.error(f"Error loading {COMPONENTS_PATH}/Components.json")
        raise Exception(f"Error loading components file")

    PATH = args.path

    #Check if PATH has 4 elements
    if len(PATH) != 4:
        logger.error("Invalid number of arguments for path. Expected: CLASS SUBCLASS TYPE SUBTYPE")
        raise Exception("Invalid number of arguments for path. Expected: CLASS SUBCLASS TYPE SUBTYPE")

    CLASS    = PATH[0]
    SUBCLASS = PATH[1]
    TYPE     = PATH[2]
    SUBTYPE  = PATH[3]

    try:
        components[CLASS][SUBCLASS]
    except KeyError:
        logger.error(f"Component {CLASS},{SUBCLASS},{TYPE},{SUBTYPE} not found")
        raise Exception(f"Component not found")

    for typ,subTyp,f in components[CLASS][SUBCLASS]:
        if typ == TYPE and subTyp == SUBTYPE:
            break

    fileFolder = UAMMD_STRUCTURED_PATH+"/"+CLASS+"/"+SUBCLASS+"/"+TYPE
    filePath   = fileFolder+"/"+f

    logger.info(f"Path of component {CLASS} {SUBCLASS} {TYPE} {SUBTYPE}: {filePath}")
