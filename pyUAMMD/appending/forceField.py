import copy

from deepdiff import DeepDiff

import logging

from ..utils.update.groups import updateComponentSetGroupsLists
from ..utils.update.groups import appendGroups

from ..utils.update.ids import updateComponentSetIds

def appendForceField(topologyKey,
                     structureKey,
                     forceFieldKey,
                     sim,sim2app,mode,idOffset,modeOffset,availGroupTypes,
                     id_labels,id_list_labels,type_labels,ignoredEntriesType=[]):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending force field
    #Force field structure:
    #   forceFieldKey: {
    #       "groups": {...} #Groups list, optional
    #       "interaction1": {
    #           "type": ["Class1","SubClass1"],
    #           "parameters": {
    #               "param1": "value",
    #               "param2": "value2",
    #               ...
    #           }
    #           "labels": ["label1","label2",...],
    #           "data": [[...],
    #                    [...],
    #                    ...]
    #       },
    #       "interaction2": {
    #           "type": ["Class2","SubClass2"],
    #           "parameters": {
    #               "param1": "value",
    #               "param2": "value2",
    #               ...
    #           }
    #           "labels": ["label1","label2",...],
    #           "data": [[...],
    #                    [...],
    #                    ...]
    #       },
    #       ...,
    #   }

    # Four cases:
    # 1) forceField NOT in sim and forceField NOT in sim2app
    # 2) forceField NOT in sim and forceField in sim2app
    # 3) forceField in sim and forceField NOT in sim2app
    # 4) forceField in sim and forceField in sim2app

    forceFieldInSim     = (topologyKey in sim)
    forceFieldInSim2app = (topologyKey in sim2app)

    if forceFieldInSim:
        forceFieldInSim     = (forceFieldKey in sim[topologyKey])
    if forceFieldInSim2app:
        forceFieldInSim2app = (forceFieldKey in sim2app[topologyKey])

    structureInSim     = (topologyKey in sim)
    if structureInSim:
        structureInSim     = (structureKey in sim[topologyKey])

    structureInSim2app = (topologyKey in sim2app)
    if structureInSim2app:
        structureInSim2app    = (structureKey in sim2app[topologyKey])

    if (forceFieldInSim or forceFieldInSim2app) and not (structureInSim or structureInSim2app):
        logger.warning("Structure must be present in at least one simulation")
        logger.error("Cannot merge forceField without any structure")
        raise Exception("Cannot merge forceField without structure")

    #Move ignored entries to ignoredEntries
    ignoredEntries = {}
    ignoredEntriesNames = []
    if forceFieldInSim2app:
        for entry in sim2app[topologyKey][forceFieldKey]:
            entryType = sim2app[topologyKey][forceFieldKey][entry]["type"][0]
            if entryType in ignoredEntriesType:
                ignoredEntriesNames.append(entry)

        for entry in ignoredEntriesNames:
            logger.debug("Moving ignored entry %s to ignoredEntries",entry)
            #Remove from sim2app (pop) and store in ignoredEntries
            ignoredEntries[entry] = sim2app[topologyKey][forceFieldKey].pop(entry)

    # 1) forceField NOT in sim and forceField NOT in sim2app
    if not forceFieldInSim and not forceFieldInSim2app:
        #Do nothing
        pass
    # 2) forceField NOT in sim and forceField in sim2app
    elif not forceFieldInSim and forceFieldInSim2app:
        #Check if the sim2app have particles defined.
        if structureInSim2app:
            #In this case, we first update the groups in sim2app and then de ids in the remainder entries
            sim[topologyKey][forceFieldKey] = updateComponentSetIds(mode,
                                                                  updateComponentSetGroupsLists(mode,
                                                                                                sim2app[topologyKey][forceFieldKey],
                                                                                                idOffset,modeOffset,availGroupTypes),
                                                                  idOffset,id_labels,id_list_labels)
        else:
            sim[topologyKey][forceFieldKey] = copy.deepcopy(sim2app[topologyKey][forceFieldKey])

    # 3) forceField in sim and forceField NOT in sim2app
    elif forceFieldInSim and not forceFieldInSim2app:
        #Do nothing
        pass
    # 4) forceField in sim and forceField in sim2app
    elif forceFieldInSim and forceFieldInSim2app:

        #Check if the sim2app have particles defined.
        if structureInSim2app:
            #We append (after updating) the groupsLists in sim2app[forceFieldKey] to sim[forceFieldKey]
            updatedForceField2app = updateComponentSetIds(mode,sim2app[topologyKey][forceFieldKey],idOffset,id_labels,id_list_labels)
            appendGroups(mode,sim[topologyKey][forceFieldKey],updatedForceField2app,idOffset,modeOffset,availGroupTypes)
        else:
            updatedForceField2app = copy.deepcopy(sim2app[topologyKey][forceFieldKey])
            appendGroups(mode,sim[topologyKey][forceFieldKey],updatedForceField2app,0,0,availGroupTypes)

        forceFieldSim     = [k for k in sim[topologyKey][forceFieldKey].keys() if sim[topologyKey][forceFieldKey][k]["type"] != ["Groups","GroupsList"]]
        forceFieldSim2app = [k for k in sim2app[topologyKey][forceFieldKey].keys() if sim2app[topologyKey][forceFieldKey][k]["type"] != ["Groups","GroupsList"]]

        for ff2app in forceFieldSim2app:
            if ff2app not in forceFieldSim:
                sim[topologyKey][forceFieldKey][ff2app] = copy.deepcopy(updatedForceField2app[ff2app])
            else:

                equalParam = not DeepDiff(sim[topologyKey][forceFieldKey][ff2app],
                                          updatedForceField2app[ff2app],
                                          ignore_order=True,report_repetition=True,
                                          exclude_paths=["root['data']","root['labels']"]) #Note data is excluded from comparison

                if equalParam:

                    #Adding force field entry to sim

                    if "data" not in sim[topologyKey][forceFieldKey][ff2app] or "data" not in updatedForceField2app[ff2app]:
                        continue

                    dataSim     = sim[topologyKey][forceFieldKey][ff2app]["data"]
                    dataSim2app = updatedForceField2app[ff2app]["data"]

                    isBatchIdDependent   = False
                    isSim2appIdDependent = False
                    for label in sim[topologyKey][forceFieldKey][ff2app]["labels"]:
                        if label in id_labels or label in id_list_labels:
                            isBatchIdDependent = True
                            break
                    for label in updatedForceField2app[ff2app]["labels"]:
                        if label in id_labels or label in id_list_labels:
                            isSim2appIdDependent = True
                            break

                    if isBatchIdDependent != isSim2appIdDependent:
                        logger.error("Force field entry \"{}\" is id dependent in one of the simulations but not in the other".format(ff2app))
                        raise Exception("ForceField entry is id dependent in one of the simulations but not in the other")
                    else:
                        isIdDependent = isBatchIdDependent

                    #Check equal data
                    if isIdDependent and (mode == "batchId"):
                        #We assume data have changed, to force that equalData = False, this avoids the need to check if the data is equal
                        valChanged = True
                    else:
                        #equalData = not DeepDiff(dataSim,dataSim2app,ignore_order=True,report_repetition=True)
                        ddiff = DeepDiff(dataSim,dataSim2app,ignore_order=True,report_repetition=True)

                        valChanged = ("values_changed" in ddiff.keys())

                    if valChanged:
                        equalData = False
                    else:
                        #Data can not be equal but dataSim2app can be a subset of dataSim or viceversa
                        added   = ddiff.get("iterable_item_added",[])
                        removed = ddiff.get("iterable_item_removed",[])

                        nAdded   = len(added)
                        nRemoved = len(removed)

                        if nAdded == 0 and nRemoved == 0:
                            equalData = True
                        if nAdded > 0 and nRemoved == 0:
                            #dataSim is a subset of dataSim2app, update dataSim
                            sim[topologyKey][forceFieldKey][ff2app]["data"] = copy.deepcopy(dataSim2app)
                            equalData = True
                        if nAdded == 0 and nRemoved > 0:
                            #dataSim2app is a subset of dataSim, do nothing
                            equalData = True
                        if nAdded > 0 and nRemoved > 0:
                            equalData = False

                    if equalData:
                        pass
                    else:

                        isSimTypeDependent     = False
                        isSim2appTypeDependent = False
                        for label in sim[topologyKey][forceFieldKey][ff2app]["labels"]:
                            if label in type_labels:
                                isSimTypeDependent = True
                                break
                        for label in updatedForceField2app[ff2app]["labels"]:
                            if label in type_labels:
                                isSim2appTypeDependent = True
                                break

                        if isSimTypeDependent != isSim2appTypeDependent:
                            logger.error("Force field entry \"{}\" is type dependent in one of the simulations but not in the other".format(ff2app))
                            raise Exception("ForceField entry is type dependent in one of the simulations but not in the other")
                        else:
                            isTypeDependent = isSimTypeDependent

                        if isIdDependent or not isTypeDependent or (mode == "modelId"):
                            sim[topologyKey][forceFieldKey][ff2app]["data"].extend(updatedForceField2app[ff2app]["data"])
                        else:
                            if(mode == "batchId"):

                                isBatchIdPresentSim     = ("batchId" in sim[topologyKey][forceFieldKey][ff2app]["labels"])
                                isBatchIdPresentSim2app = ("batchId" in updatedForceField2app[ff2app]["labels"])

                                if not isBatchIdPresentSim and not isBatchIdPresentSim2app:
                                    sim[topologyKey][forceFieldKey][ff2app]["labels"].append("batchId")
                                    updatedForceField2app[ff2app]["labels"].append("batchId")

                                    for i in range(len(sim[topologyKey][forceFieldKey][ff2app]["data"])):
                                        sim[topologyKey][forceFieldKey][ff2app]["data"][i].append(0)
                                    for i in range(len(updatedForceField2app[ff2app]["data"])):
                                        updatedForceField2app[ff2app]["data"][i].append(1)

                                    sim[topologyKey][forceFieldKey][ff2app]["data"].extend(updatedForceField2app[ff2app]["data"])

                                elif isBatchIdPresentSim and not isBatchIdPresentSim2app:
                                    updatedForceField2app[ff2app]["labels"].append("batchId")

                                    for i in range(len(updatedForceField2app[ff2app]["data"])):
                                        updatedForceField2app[ff2app]["data"][i].append(modeOffset)

                                    sim[topologyKey][forceFieldKey][ff2app]["data"].extend(updatedForceField2app[ff2app]["data"])

                                elif not isBatchIdPresentSim and isBatchIdPresentSim2app:
                                    logger.error("BatchId not present in reference simulation, but present in simulation to append.\
                                                 Force field entry: \"{}\"".format(ff2app))
                                    raise Exception("BatchId not present in reference simulation, but present in simulation to append")
                                else:
                                    #Both are present
                                    logger.error("BatchId present in both simulations. Force field entry: \"{}\"".format(ff2app))
                                    raise Exception("BatchId present in both simulations")


                            else:
                                logger.error(f"This should not happen, please report it to the developers. \
                                             isIdDependent: {isIdDependent}, isTypeDependent: {isTypeDependent}, mode: {mode}")
                                raise Exception("This should not happen, please report it to the developers")

                else:
                    logger.error("Only force field entries which differ in data can be appended. Force field entry: \"{}\"".format(ff2app))
                    raise Exception("Only force field entries which differ in data can be appended")

    #If ignoreEntries is not empty, move back entries to sim2app
    if len(ignoredEntriesNames) > 0:
        for entry in ignoredEntriesNames:
            logger.debug("Moving back entry {} to simulation to append".format(entry))
            sim2app[topologyKey][forceFieldKey][entry] = ignoredEntries.pop(entry)

    #Force field appended
    #########################
