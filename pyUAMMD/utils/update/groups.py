import copy
import logging

from deepdiff import DeepDiff

def getGroupsListEntriesFromComponentSet(componentSet):
    """
    Returns a list with all the entries which class is "Groups" and subClass "GroupsList"
    """
    groups = []
    for name,entry in componentSet.items():
        if "type" in entry:
            if entry["type"][0] == "Groups":
                if entry["type"][1] == "GroupsList":
                    groups.append(name)

    return groups

def updateComponentSetGroupsLists(mode,componentSet,idOffset,modeOffset,availGroupTypes):
    """
    This function take a component set. It looks for all the groups list entries
    in the component set. Then, if some groups list is found, it updates the
    group definition according to the mode and values of idOffset and modeOffset.

    Returns a copy of the component set with the updated groups.
    """

    logger = logging.getLogger("pyUAMMD")

    updatedComponentSet = copy.deepcopy(componentSet)

    groupsLists = getGroupsListEntriesFromComponentSet(updatedComponentSet)

    #If not, doing nothing
    if len(groupsLists) == 0:
        return updatedComponentSet

    for name in groupsLists:
        logger.debug("Updating groups list \"{}\"".format(name))

        groupDeclarations = updatedComponentSet[name]["data"]

        for groupName,groupType,selection in groupDeclarations:
            if groupType not in availGroupTypes:
                logger.error("Group type \"{}\" is not supported.".format(groupType))
                raise Exception("Group type not supported.")
            else:
                if groupType in ["Ids","notIds"]:
                    for n,_ in enumerate(selection):
                        selection[n] += idOffset
                elif groupType in ["Types","notTypes"]:
                    #Do nothing
                    pass
                elif groupType in ["ModelIds","notModelIds"]:
                    if mode == "modelId":
                        for n,_ in enumerate(selection):
                            selection[n] += modeOffset
                    else:
                        #Do nothing
                        pass
                elif groupType in ["BatchIds","notBatchIds"]:
                    if mode == "batchId":
                        for n,_ in enumerate(selection):
                            selection[n] += modeOffset
                    else:
                        #Do nothing
                        pass
                elif groupType in ["ModelIdsBatchIds","notModelIdsBatchIds"]:
                    if mode == "modelId":
                        for n,_ in enumerate(selection):
                            selection[n][0] += modeOffset
                    elif mode == "batchId":
                        for n,_ in enumerate(selection):
                            selection[n][1] += modeOffset
                    else:
                        #Do nothing
                        pass

    return updatedComponentSet


def appendGroups(mode,componentRefSet,component2AppendSet,idOffset,modeOffset,availGroupTypes):
    """
    This function takes two component set. It updates all the groupsList entries found
    in the component set to append according to the mode and values of idOffset and modeOffset.
    Then it appends the groupsList entries of the component set to append to the component set
    of reference.

    It does not return anything. Reference component set is updated.

    Rules:
    *No groups list in reference component set:
    -- No groups list in append component set -> Do nothing
    -- Groups list in append component set -> Update all components of the group list type in the set of components to be added and append them to reference
    *Groups list in reference component set:
    -- No groups list in append component set -> Do nothing
    -- Groups list in append component set ->
            1) Update all components of the group list type in the set of components to be added
            2) If the name of some group list type in the set of components to be added
               match some name of some group list type in the reference set of components ->
               1) If both group list types are the same and selection is the same -> Do nothing
               2) If both group list types are the same and selection is different -> Try to append the group list type to the reference set of components
               3) If both group list types are different -> Error

    """

    logger = logging.getLogger("pyUAMMD")

    #Check if a entry of class "Groups" and subClass "GroupsList" is present
    groupsListsRef     = getGroupsListEntriesFromComponentSet(componentRefSet)
    groupsLists2Append = getGroupsListEntriesFromComponentSet(component2AppendSet)

    #If not, doing nothing
    if len(groupsListsRef) == 0:
        if len(groupsLists2Append) == 0:
            return
        else:
            updatedComponent2AppendSet = updateComponentSetGroupsLists(mode,component2AppendSet,idOffset,modeOffset,availGroupTypes)
            for name in groupsLists2Append:
                componentRefSet[name] = copy.deepcopy(updatedComponent2AppendSet[name])
    else:
        if len(groupsLists2Append) == 0:
            return
        else:
            updatedComponent2AppendSet = updateComponentSetGroupsLists(mode,component2AppendSet,idOffset,modeOffset,availGroupTypes)

            for name in groupsLists2Append:
                if name in groupsListsRef:

                    refGroupDeclarationNames      = [x[0] for x in componentRefSet[name]["data"]]
                    _2AppendGroupDeclarationNames = [x[0] for x in component2AppendSet[name]["data"]]

                    #Check groups names are unique in each set
                    if len(list(set(refGroupDeclarationNames))) != len(refGroupDeclarationNames):
                        logger.error("There are repeated group declaration names in the component set of reference.")
                        raise Exception("Repeated group declaration names.")
                    if len(list(set(_2AppendGroupDeclarationNames))) != len(_2AppendGroupDeclarationNames):
                        logger.error("There are repeated group declaration names in the component set to append.")
                        raise Exception("Repeated group declaration names.")

                    #Check there are no repeteaded names
                    if len(refGroupDeclarationNames) != len(set(refGroupDeclarationNames)):
                        logger.error("There are repeated group declaration names in the component set of reference.")
                        raise Exception("Repeated group declaration names.")
                    if len(_2AppendGroupDeclarationNames) != len(set(_2AppendGroupDeclarationNames)):
                        logger.error("There are repeated group declaration names in the component set to append.")
                        raise Exception("Repeated group declaration names.")

                    for gName2Append in _2AppendGroupDeclarationNames:
                        if gName2Append in refGroupDeclarationNames:
                            gRef = [gDecl for gDecl in componentRefSet[name]["data"] if gDecl[0] == gName2Append]
                            gApp = [gDecl for gDecl in updatedComponent2AppendSet[name]["data"] if gDecl[0] == gName2Append]

                            # The previous list must have only one element since the names are unique
                            gRef = gRef[0]
                            gApp = gApp[0]

                            equal = not DeepDiff(gRef,gApp,ignore_order=True)
                            if equal:
                                continue
                            else:
                                typeRef = gRef[1]
                                typeApp = gApp[1]
                                if typeRef != typeApp:
                                    logger.error(f"Group types are different for group {gName2Append}.")
                                    raise Exception("Group types are different.")
                                else:
                                    #Merge the two lists
                                    #Find group index in "data" list
                                    idx = [i for i,x in enumerate(componentRefSet[name]["data"]) if x[0] == gName2Append][0]
                                    #Append the new selection
                                    componentRefSet[name]["data"][idx][2] += gApp[2]
                                    #Remove duplicates
                                    componentRefSet[name]["data"][idx][2] = list(set(componentRefSet[name]["data"][idx][2]))
                        else:
                            g = [gDecl for gDecl in updatedComponent2AppendSet[name]["data"] if gDecl[0] == gName2Append]

                            componentRefSet[name]["data"].append(g[0])
                else:
                    for name in groupsLists2Append:
                        componentRefSet[name] = copy.deepcopy(updatedComponent2AppendSet[name])
