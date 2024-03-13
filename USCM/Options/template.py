import sys
import os

import logging

import json

def template(args,
             COMPONENTS_PATH,TEMPLATES_PATH):

    logger = logging.getLogger("USCM")

    logger.info(f"Loading components file {COMPONENTS_PATH}/Components.json")
    try:
        with open(COMPONENTS_PATH+"/Components.json") as f:
            components = json.load(f)
    except:
        logger.error(f"Error loading {COMPONENTS_PATH}/Components.json")
        raise Exception(f"Error loading components file")

    #Try to load the Components.json file
    logger.info(f"Loading templates file {TEMPLATES_PATH}/Templates.json\n")

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

    logger.debug(f"Checking templates consistency")

    #Check all entries in Templates.json are also in Components.json
    for CLASS in templates:

        try:
            components[CLASS]
        except KeyError:
            logger.error(f"Class \"{CLASS}\" is in Templates.json but not in Components.json")
            raise Exception(f"Class in Templates.json but not in Components.json")

        for SUB_CLASS in templates[CLASS]:

            try:
                components[CLASS][SUB_CLASS]
            except KeyError:
                logger.error(f"Class,SubClass \"{CLASS},{SUB_CLASS}\" is in Templates.json but not in Components.json")
                raise Exception(f"Class,SubClass in Templates.json but not in Components.json")

            for TYPE in templates[CLASS][SUB_CLASS]:
                found = False
                for typ,subTyp,f in components[CLASS][SUB_CLASS]:
                    if typ == TYPE:
                        found = True

                if not found:
                    logger.error(f"Class,SubClass,Type \"{CLASS},{SUB_CLASS},{TYPE}\" is in Templates.json but not in Components.json")
                    raise Exception(f"Class,SubClass,Type in Templates.json but not in Components.json")

    logger.info("Creating a new component template...")

    if len(args.template) != 4:
        logger.error("Invalid number of arguments, please provide: CLASS, SUBCLASS, TYPE, SUBTYPE")
        raise Exception("Invalid number of arguments")

    CLASS    = args.template[0]
    SUBCLASS = args.template[1]
    TYPE     = args.template[2]
    SUBTYPE  = args.template[3]

    if args.variant:
        VARIANT = args.variant[0]
    else:
        VARIANT = None

    requestedTemplate = [CLASS,SUBCLASS,TYPE,VARIANT]
    #Found closest template [CLASS,SUBCLASS,TYPE,VARIANT] in the list of available templates
    points = [0 for i in range(len(availTemplates))]
    for i,template in enumerate(availTemplates):
        for j,level in enumerate(template):
            if level == requestedTemplate[j]:
                points[i] += 1

    #print points per available template
    for i,template in enumerate(availTemplates):
        logger.debug(f"Template {template} has {points[i]} points")

    #Selected template is the template with the highest number of points
    #If two templates have the same number of points, we select the shortest one
    maxPoints = max(points)
    selectedTemplate = None
    for i,template in enumerate(availTemplates):
        if points[i] == maxPoints:
            if selectedTemplate is None:
                selectedTemplate = template
            else:
                if len(template) < len(selectedTemplate):
                    selectedTemplate = template

    logger.info(f"Selected template: {selectedTemplate}")

    #Check if the template file exists
    templateFilePath = os.path.join(TEMPLATES_PATH,"TemplatesCode/")
    for level in selectedTemplate:
        templateFilePath += f"{level}_"
    templateFilePath = templateFilePath[:-1] + ".cuh"
    if not os.path.exists(templateFilePath):
        logger.error(f"Template file {templateFilePath} not found")
        raise Exception(f"Template file not found")
    else:
        with open(templateFilePath) as f:
            template = f.read()
            template = template.replace("__CLASS__",CLASS)
            template = template.replace("__SUBCLASS__",SUBCLASS)
            template = template.replace("___TYPE___","_"+TYPE.upper()+"_")
            template = template.replace("__TYPE__",TYPE)

            template = template.replace("___SUBTYPE___","_"+SUBTYPE.upper()+"_")
            template = template.replace("__SUBTYPE__",SUBTYPE)

            #Save the file in the current directory
            #Check if file SUBTYPE.cuh exists
            if os.path.exists(f"{SUBTYPE}.cuh"):
                logger.error(f"File {SUBTYPE}.cuh already exists")
                raise Exception(f"File {SUBTYPE}.cuh already exists")
            else:
                with open(f"{SUBTYPE}.cuh", "w") as f:
                    f.write(template)
                    logger.info(f"Component {SUBTYPE} created")
