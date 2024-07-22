import logging

from ..common import getDataEntryLabelIndex

import copy

def updateDataEntryIds(entry,idOffset,id_labels,id_list_labels):
    """
    Expected a data entry
    Data is updated if it contains id labels.
    Data is considered id if its label is in id_labels
    Data is considered id_list if its label is in id_list_labels
    """

    logger = logging.getLogger("pyUAMMD")

    if "labels" not in entry.keys():
        return
    if "data" not in entry.keys():
        return

    indexToUpdate = []

    for l in entry["labels"]:
        if l in id_labels:
            indexToUpdate.append(getDataEntryLabelIndex(entry,l))

    for d in entry["data"]:
        for i in indexToUpdate:
            d[i] += idOffset

    indexToUpdate = []

    for l in entry["labels"]:
        if l in id_list_labels:
            indexToUpdate.append(getDataEntryLabelIndex(entry,l))

    for d in entry["data"]:
        for i in indexToUpdate:
            for n,_ in enumerate(d[i]):
                d[i][n] += idOffset

def updateComponentSetIds(mode,componentSet,idOffset,id_labels,id_list_labels):
    """
    This function take a component set. It looks for all the entries that
    contains data refered to id. Then, if some id is found, it updates the
    id according to the mode and values of idOffset.

    Returns a copy of the component set with the updated ids.
    """

    logger = logging.getLogger("pyUAMMD")

    updatedComponentSet = copy.deepcopy(componentSet)

    for dataEntry in updatedComponentSet.values():
        #If entry is not a data entry, updateDataEntryIds will do nothing
        updateDataEntryIds(dataEntry,idOffset,id_labels,id_list_labels)

    return updatedComponentSet
