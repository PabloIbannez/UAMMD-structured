import os
import sys

import logging

import json

from colorama import init as colorama_init
from colorama import Fore
from colorama import Style

def lst(args,
        COMPONENTS_PATH,DATA_PATH,TEMPLATES_PATH):

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

    if args.list[0] == "components":
        logger.info("Listing all the components available in UAMMD Structured...")
        for component in components:
            #Each component is a dictionary, we can iterate over the keys
            logger.info(f"{Fore.GREEN}Component class: {component}")
            for key in components[component]:
                componentList = components[component][key]
                #Print every element in the list in a new line
                logger.info(f"{Fore.BLUE}  Component subClass: {key}")

                typs = set([element[0] for element in componentList])
                for t in typs:
                    logger.info(f"    * {Fore.MAGENTA}Type: {t}")
                    for element in componentList:
                        if element[0] == t:
                            logger.info(f"      - SubType: {element[1]}")

    elif args.list[0] == "data":
        logger.info("Listing all the data available in UAMMD Structured...")
        for d in data:
            #Each component is a dictionary, we can iterate over the keys
            logger.info(f"{Fore.GREEN}Data: {d}")
            for dlist in data[d]:

                name   = dlist[0]
                getset = dlist[1]
                typ    = dlist[2]

                try:
                    default = dlist[3]
                except IndexError:
                    default = None

                logger.info(f"{Fore.BLUE}  Name: {name} {Fore.RESET}")
                logger.info(f"    - Getter/Setter: {getset}")
                logger.info(f"    - Type: {typ}")
                logger.info(f"    - Default: {default}")


    elif args.list[0] == "templates":
        try:
            with open(TEMPLATES_PATH+"/Templates.json") as f:
                templates = json.load(f)

            availTemplates = []
            for CLASS in templates:
                availTemplates.append([CLASS])
                for SUBCLASS in templates[CLASS]:
                    availTemplates.append([CLASS,SUBCLASS])
                    for TYPE in templates[CLASS][SUBCLASS]:
                        availTemplates.append([CLASS,SUBCLASS,TYPE])

                        variants = templates[CLASS][SUBCLASS][TYPE]
                        for variant in variants:
                            availTemplates.append([CLASS,SUBCLASS,TYPE,variant])

        except FileNotFoundError:
            logger.error(f"Error loading {TEMPLATES_PATH}/Templates.json")
            raise Exception(f"Error loading templates file")

        logger.info("Listing all the templates available in UAMMD Structured...")
        for template in templates:
            logger.info(f"{Fore.GREEN}Template class: {template}")
            for key in templates[template]:
                templateList = templates[template][key]
                logger.info(f"{Fore.BLUE}  Template subClass: {key}")
                for element in templateList:
                    logger.info(f"    * {Fore.MAGENTA}type: {element}")
                    for variant in templateList[element]:
                        logger.info(f"      - Variant: {variant}")

        #Check if all templates are present
        logger.info(f"")
        logger.info(f"Checking if all templates files are present...")

        for template in availTemplates:
            templateFilePath = os.path.join(TEMPLATES_PATH,"TemplatesCode/")
            for level in template:
                templateFilePath += f"{level}_"
            templateFilePath = templateFilePath[:-1] + ".cuh"
            if not os.path.exists(templateFilePath):
                logger.warning(f"Template file {templateFilePath} {Fore.RED}NOT FOUND{Fore.RESET}")

    else:
        logger.error("Invalid argument for list, available options are: components, templates, data")
        raise Exception("Invalid argument for list, available options are: components, templates, data")

