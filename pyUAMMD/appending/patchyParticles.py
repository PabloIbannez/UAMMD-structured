import copy

from deepdiff import DeepDiff

import logging

from ..utils.common import getDataEntryLabelIndex

from ..utils.update.groups import updateComponentSetGroupsLists
from ..utils.update.groups import appendGroups

from ..utils.update.ids import updateComponentSetIds

from .glb import appendGlobal
from .state import appendState
from .forceField import appendForceField

def appendPatchyParticlesStructure(topologyKey,
                                   structureKey,
                                   patchesSim,patchesSim2app,mode,
                                   parentIdOffset):

    logger = logging.getLogger("pyUAMMD")

    idOffset   = 0
    modeOffset = 0

    #########################
    #Appending structure
    #Structure structure:
    #   structureKey: {
    #       "labels": ["id","type","resId","chainId","modelId","batchId"]
    #       "data": [[...],
    #                [...],
    #                ...]

    # Four cases:
    # 1) structure NOT in patchesSim and structure NOT in patchesSim2app
    # 2) structure NOT in patchesSim and structure in patchesSim2app
    # 3) structure in patchesSim and structure NOT in patchesSim2app
    # 4) structure in patchesSim and structure in patchesSim2app

    structureInPatchesSim     = (topologyKey in patchesSim)
    structureInPatchesSim2app = (topologyKey in patchesSim2app)

    if structureInPatchesSim:
        structureInPatchesSim     = (structureKey in patchesSim[topologyKey])
    if structureInPatchesSim2app:
        structureInPatchesSim2app = (structureKey in patchesSim2app[topologyKey])

    #Create offsets if reference has structure
    if structureInPatchesSim:

        #Computing id offset
        idIndex = getDataEntryLabelIndex(patchesSim[topologyKey][structureKey],"id")

        idOffset = 0
        for d in patchesSim[topologyKey][structureKey]["data"]:
            if d[idIndex] > idOffset:
                idOffset = d[idIndex]
        idOffset += 1

        modeOffset = 0
        if mode == "modelId":
            modeOffset = parentIdOffset
        elif mode == "batchId":
            if "batchId" not in patchesSim[topologyKey][structureKey]["labels"]:
                patchesSim[topologyKey][structureKey]["labels"].append("batchId")
                for d in patchesSim[topologyKey][structureKey]["data"]:
                    d.append(0)

            modeIndex = getDataEntryLabelIndex(patchesSim[topologyKey][structureKey],"batchId")
            for d in patchesSim[topologyKey][structureKey]["data"]:
                if d[modeIndex] > modeOffset:
                    modeOffset = d[modeIndex]
            modeOffset += 1
        else:
            logger.error(f"Appending mode {mode} not implemented")
            raise Exception("Appending mode not implemented")

    # 1) structure NOT in patchesSim and structure NOT in patchesSim2app
    if not structureInPatchesSim and not structureInPatchesSim2app:
        #Do nothing
        pass
    # 2) structure NOT in patchesSim and structure in patchesSim2app
    elif not structureInPatchesSim and structureInPatchesSim2app:
        if topologyKey not in patchesSim:
            patchesSim[topologyKey] = {}
        patchesSim[topologyKey][structureKey] = copy.deepcopy(patchesSim2app[topologyKey][structureKey])

        #Data merging
        for dataIndex in range(len(patchesSim[topologyKey][structureKey]["data"])):
            for label in patchesSim[topologyKey][structureKey]["labels"]:
                if label == "parentId":
                    parentIdIndex = getDataEntryLabelIndex(patchesSim[topologyKey][structureKey],"parentId")
                    patchesSim[topologyKey][structureKey]["data"][dataIndex][parentIdIndex] += parentIdOffset

        idOffset   = 0
        if mode == "modelId":
            modeOffset = parentIdOffset
        elif mode == "batchId":
            modeOffset = 0
        else:
            logger.error(f"Appending mode {mode} not implemented")
            raise Exception("Appending mode not implemented")

    # 3) structure in patchesSim and structure NOT in patchesSim2app
    elif structureInPatchesSim and not structureInPatchesSim2app:
        #Do nothing
        pass
    # 4) structure in patchesSim and structure in patchesSim2app
    elif structureInPatchesSim and structureInPatchesSim2app:

        #Handle labels
        patchesSimLabels     = patchesSim[topologyKey][structureKey]["labels"]
        patchesSim2appLabels = patchesSim2app[topologyKey][structureKey]["labels"]

        availLabels = ["id","type","parentId","batchId"]

        for label in patchesSimLabels:
            if label not in availLabels:
                logger.error(f"Label {label} not available for structure.")
                raise Exception("Label not available for structure.")

        for label in patchesSim2appLabels:
            if label not in availLabels:
                logger.error(f"Label {label} not available for structure.")
                raise Exception("Label not available for structure.")

        for label in availLabels:

            labelInPatchesSim     = (label in patchesSimLabels)
            labelInPatchesSim2app = (label in patchesSim2appLabels)

            if not labelInPatchesSim and not labelInPatchesSim2app:
                #Do nothing
                pass
            elif not labelInPatchesSim and labelInPatchesSim2app:
                patchesSim[topologyKey][structureKey]["labels"].append(label)
                for d in patchesSim[topologyKey][structureKey]["data"]:
                    d.append(0)
            elif labelInPatchesSim and not labelInPatchesSim2app:
                patchesSim2app[topologyKey][structureKey]["labels"].append(label)
                for d in patchesSim2app[topologyKey][structureKey]["data"]:
                    d.append(0)
            elif labelInPatchesSim and labelInPatchesSim2app:
                #Do nothing
                pass

        #Data merging
        for structInfo in patchesSim2app[topologyKey][structureKey]["data"]:
            d = []
            for label in patchesSim[topologyKey][structureKey]["labels"]:
                if label == "id":
                    d.append(structInfo[patchesSim2appLabels.index(label)] + idOffset)
                elif label == "batchId" and mode == "batchId":
                    d.append(structInfo[patchesSim2appLabels.index(label)] + modeOffset)
                elif label == "parentId":
                    d.append(structInfo[patchesSim2appLabels.index(label)] + parentIdOffset)
                else:
                    d.append(structInfo[patchesSim2appLabels.index(label)])
            patchesSim[topologyKey][structureKey]["data"].append(d)

    #Structure appended
    #########################

    return idOffset,modeOffset

def appendPatchyParticles(topologyKey,
                          structureKey,
                          forceFieldKey,
                          sim,sim2app,mode,
                          availGroupTypes,
                          id_labels,id_list_labels,type_labels):

    logger = logging.getLogger("pyUAMMD")

    forceFieldInSim     = (topologyKey in sim)
    forceFieldInSim2app = (topologyKey in sim2app)

    if forceFieldInSim:
        forceFieldInSim     = (forceFieldKey in sim[topologyKey])
    if forceFieldInSim2app:
        forceFieldInSim2app = (forceFieldKey in sim2app[topologyKey])

    if not forceFieldInSim and not forceFieldInSim2app:
        return

    patchyParticlesInSim     = []
    patchyParticlesInSim2app = []

    if forceFieldInSim:
        for ffEntry in sim[topologyKey][forceFieldKey]:
            if sim[topologyKey][forceFieldKey][ffEntry]["type"][0] == "PatchyParticles":
                patchyParticlesInSim.append(ffEntry)

    if forceFieldInSim2app:
        for ffEntry in sim2app[topologyKey][forceFieldKey]:
            if sim2app[topologyKey][forceFieldKey][ffEntry]["type"][0] == "PatchyParticles":
                patchyParticlesInSim2app.append(ffEntry)

    patchyParticlesEntries = []
    for ppEntry in patchyParticlesInSim:
        if ppEntry in patchyParticlesInSim2app:
            patchyParticlesEntries.append([ppEntry,ppEntry])
            patchyParticlesInSim2app.remove(ppEntry)
        else:
            patchyParticlesEntries.append([ppEntry,None])

    for ppEntry in patchyParticlesInSim2app:
        patchyParticlesEntries.append([None,ppEntry])

    for patchyParticlesInSimName,patchyParticlesInSim2appName in patchyParticlesEntries:

        patchyParticlesInSim = (patchyParticlesInSimName is not None)
        patchyParticlesInSim2app = (patchyParticlesInSim2appName is not None)

        # Four cases:
        # 1) patchyParticles NOT in sim and patchyParticles NOT in sim2app
        # 2) patchyParticles NOT in sim and patchyParticles in sim2app
        # 3) patchyParticles in sim and patchyParticles NOT in sim2app
        # 4) patchyParticles in sim and patchyParticles in sim2app

        # 1) patchyParticles NOT in sim and patchyParticles NOT in sim2app
        if not patchyParticlesInSim and not patchyParticlesInSim2app:
            #Do nothing
            pass
        # 2) patchyParticles NOT in sim and patchyParticles in sim2app
        elif not patchyParticlesInSim and patchyParticlesInSim2app:

            patchyParticlesInSimName = patchyParticlesInSim2appName

            if topologyKey not in sim:
                sim[topologyKey] = {}

            if forceFieldKey not in sim[topologyKey]:
                sim[topologyKey][forceFieldKey] = {}

            sim[topologyKey][forceFieldKey][patchyParticlesInSimName] = {}

            sim[topologyKey][forceFieldKey][patchyParticlesInSimName]["type"] = \
            copy.deepcopy(sim2app[topologyKey][forceFieldKey][patchyParticlesInSimName]["type"])

            patchyParticlesInSim = True # This ensures that the last case is executed

        # 3) patchyParticles in sim and patchyParticles NOT in sim2app
        elif patchyParticlesInSim and not patchyParticlesInSim2app:
            #Do nothing
            pass


        # 4) patchyParticles in sim and patchyParticles in sim2app
        if patchyParticlesInSim and patchyParticlesInSim2app: # Note not elif !!!!

            patchesSim     = sim[topologyKey][forceFieldKey][patchyParticlesInSim2appName]
            patchesSim2app = sim2app[topologyKey][forceFieldKey][patchyParticlesInSim2appName]

            appendGlobal("patchesGlobal",patchesSim,patchesSim2app,mode)

            #Structure requires a special treatment

            patchesIdOffset,patchesModeOffset = appendPatchyParticlesStructure("patchesTopology",
                                                                               "structure",
                                                                               patchesSim,patchesSim2app,mode,
                                                                               sim.getNumberOfParticles())

            appendState("patchesTopology","structure",
                        "patchesState",patchesSim,patchesSim2app,mode,patchesIdOffset,patchesModeOffset)

            appendForceField("patchesTopology","structure",
                             "forceField",
                             patchesSim,patchesSim2app,mode,patchesIdOffset,patchesModeOffset,availGroupTypes,
                             id_labels, id_list_labels,type_labels,
                             ignoredEntriesType=["PatchyParticles"])
