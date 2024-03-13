import os
import logging
import shutil
import json
import jsbeautifier

def copyFilesInFolder(originPath, destinationPath, logger):
    filesOrigin  = os.listdir(originPath)
    filesDestination = os.listdir(destinationPath)
    for file_origin in filesOrigin:
        if file_origin == "Components.json":
            continue
        pathToFileOrigin = os.path.join(originPath, file_origin)
        if os.path.isdir(pathToFileOrigin): #If the folder exists copy recursively the directories inside it, else copy the full folder
            pathToFolderDestination = os.path.join(destinationPath, file_origin)
            if os.path.exists(pathToFolderDestination):
                copyFilesInFolder(pathToFileOrigin, pathToFolderDestination, logger)
            else:
                pathFolderInDestination = os.path.join(destinationPath, file_origin)
                shutil.copytree(pathToFileOrigin, pathFolderInDestination)
        else: #Copy the file to the destinationPath, only if the file do not exists
            if file_origin in filesDestination:
                logger.warning(f"The file {file_origin} already exists, it will not be copied.")
            else:
                pathToFileDestination = os.path.join(destinationPath, file_origin)
                shutil.copy2(pathToFileOrigin, pathToFileDestination)


def combineTwoDictionaries(dict1, dict2, logger):
    # Merge components_extension into components_uammd
    for key, value in dict2.items():
        # If key exists in components_1 and both are lists, extend the list
        if key in dict1:
            if type(dict1[key]) is not type(value):
                logger.error(f"The variable type corresponding to the key {key} mistmach.")
                raise Exception(f"Non compatible variable type.")
            if isinstance(dict1[key], list) and isinstance(value, list):
                dict1[key].extend(value)
            elif isinstance(dict1[key], dict):
                dict1[key] = combineTwoDictionaries(dict1[key], value, logger)
            else:
                logger.error(f"Unexpected variable type if {key}")
                raise Exception(f"Unexpected varible type.")
        else:
            dict1[key] = value
    return dict1


def mergeComponentsFiles(COMPONENTS_PATH, EXTENSION_PATH, logger):
    # Load the JSON files
    with open(f'{COMPONENTS_PATH}/Components.json', 'r') as f:
        components_uammd = json.load(f)

    with open(f'{EXTENSION_PATH}/Components.json', 'r') as f:
        components_extension = json.load(f)

    components_uammd_new = combineTwoDictionaries(components_uammd, components_extension, logger)
    with open(f'{COMPONENTS_PATH}/Components.json', 'w') as f:
        f.write(jsbeautifier.beautify(json.dumps(components_uammd_new)))


def extensions(args,
               COMPONENTS_PATH,
               UAMMD_STRUCTURED_PATH):

    logger = logging.getLogger("USCM")

    logger.info("Adding an extension ...")

    EXTENSION_FOLDER = args.extensions[0]
    #Check if the extension path is a git repository
    #If so create a folder called TMP and clone the repository there
    #Then do EXTENSION_FOLDER = TMP_EXTENSION_FOLDER

    #Check if path exits
    if not os.path.exists(EXTENSION_FOLDER):
        logger.error(f"Extension folder does not exist. {EXTENSION_FOLDER}")
        raise Exception(f"Extension path does not exists.")

    #If path exits check if there is a file with the name Components.json
    components_path = os.path.join(EXTENSION_FOLDER, "Components.json")
    if not os.path.exists(components_path):
        logger.error("Components.json file not found in the extension path.")
        raise Exception("Components.json file not found.")

    #Copy all the folders subfolders and files to the main folder that is in UAMMD_STRUCTURED_PATH
    mergeComponentsFiles(COMPONENTS_PATH, EXTENSION_FOLDER, logger)
    copyFilesInFolder(EXTENSION_FOLDER, UAMMD_STRUCTURED_PATH, logger)

    if EXTENSION_FOLDER == "TMP_EXTENSION_FOLDER":
        shutil.rmtree(EXTENSION_FOLDER)
    logger.info("Extension added successfully.")



