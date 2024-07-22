import logging

def isDataEntry(entry):
    """
    Check if entry is a data entry
    """
    dataEntryKeys = ["type","data","labels","parameters"]
    for key in entry.keys():
        if key not in dataEntryKeys:
            return False
    return True

def getDataEntryLabelIndex(dataEntry,label):
    """
    Expected a data entry
    It returns the index of the label
    """

    logger = logging.getLogger("pyUAMMD")

    for i,lbl in enumerate(dataEntry["labels"]):
        if lbl == label:
            return i

    logger.error(f"Required {label} but is not present.")
    raise Exception(f"Required label is not present.")


