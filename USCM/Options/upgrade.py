import os
import logging

def upgrade(args):

    logger = logging.getLogger("USCM")

    logger.info("Upgrading UAMMD-structured ...")

    FOLDER_TO_UPGRADE = args.upgrade[0]

    #Check if folder exists
    if not os.path.isdir(FOLDER_TO_UPGRADE):
        logger.error(f"Folder {FOLDER_TO_UPGRADE} does not exist.")
        raise Exception(f"Folder does not exist.")

    try:
        UAMMD_PATH            = os.environ["UAMMD_PATH"]
    except:
        logger.error("UAMMD_PATH not defined.")
        raise Exception("UAMMD_PATH not defined.")

    try:
        UAMMD_STRUCTURED_PATH = os.environ["UAMMD_STRUCTURED_PATH"]
    except:
        logger.error("UAMMD_STRUCTURED_PATH not defined.")
        raise Exception("UAMMD_STRUCTURED_PATH not defined.")

    #Check if "preamble.h" in "../UAMMD_STRUCTURED_PATH"

    import subprocess

    if not os.path.isfile(UAMMD_PATH+"/extensions/preamble.h"):
        logger.warning(f"File {UAMMD_PATH}/extensions/preamble.h does not exist.")
    else:
        if not os.path.isfile(FOLDER_TO_UPGRADE+"/preamble.h"):
            logger.warning(f"File {FOLDER_TO_UPGRADE}/preamble.h does not exist.")
        else:
            #Use rsync to update preamble.h
            logger.info("Updating preamble.h ...")
            subprocess.run(["rsync","-au","--delete",UAMMD_PATH+"/extensions/preamble.h",FOLDER_TO_UPGRADE+"/preamble.h"])

    #Check if "structured" folder exists in FOLDER_TO_UPGRADE

    if not os.path.isdir(FOLDER_TO_UPGRADE+"/structured"):
        logger.warning(f"Folder {FOLDER_TO_UPGRADE}/structured does not exist.")
    else:
        #Use rsync to update structured folder
        logger.info("Updating structured folder ...")
        subprocess.run(["rsync","-au","--delete",UAMMD_STRUCTURED_PATH+"/",FOLDER_TO_UPGRADE+"/structured/"])

    #Check if "USCM" folder exists in FOLDER_TO_UPGRADE

    if not os.path.isdir(FOLDER_TO_UPGRADE+"/USCM"):
        logger.warning(f"Folder {FOLDER_TO_UPGRADE}/USCM does not exist.")
    else:
        #Use rsync to update USCM folder
        logger.info("Updating USCM folder ...")
        subprocess.run(["rsync","-au","--delete",UAMMD_PATH+"/USCM/",FOLDER_TO_UPGRADE+"/USCM/"])
