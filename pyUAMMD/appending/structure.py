import copy

import logging

from ..utils.common import getDataEntryLabelIndex

def appendStructure(topologyKey,
                    structureKey,
                    sim,sim2app,mode):

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
    # 1) structure NOT in sim and structure NOT in sim2app
    # 2) structure NOT in sim and structure in sim2app
    # 3) structure in sim and structure NOT in sim2app
    # 4) structure in sim and structure in sim2app

    structureInSim     = (topologyKey in sim)
    structureInSim2app = (topologyKey in sim2app)

    if structureInSim:
        structureInSim     = (structureKey in sim[topologyKey])
    if structureInSim2app:
        structureInSim2app = (structureKey in sim2app[topologyKey])

    #Create offsets if reference has structure
    if structureInSim:

        #Computing id offset
        idIndex = getDataEntryLabelIndex(sim[topologyKey][structureKey],"id")

        idOffset = 0
        for d in sim[topologyKey][structureKey]["data"]:
            if d[idIndex] > idOffset:
                idOffset = d[idIndex]
        idOffset += 1

        modeOffset = 0
        if mode == "modelId":
            if "modelId" not in sim[topologyKey][structureKey]["labels"]:
                sim[topologyKey][structureKey]["labels"].append("modelId")
                for d in sim[topologyKey][structureKey]["data"]:
                    d.append(0)

            modeIndex = getDataEntryLabelIndex(sim[topologyKey][structureKey],"modelId")
            for d in sim[topologyKey][structureKey]["data"]:
                if d[modeIndex] > modeOffset:
                    modeOffset = d[modeIndex]
            modeOffset += 1
        elif mode == "batchId":
            if "batchId" not in sim[topologyKey][structureKey]["labels"]:
                sim[topologyKey][structureKey]["labels"].append("batchId")
                for d in sim[topologyKey][structureKey]["data"]:
                    d.append(0)

            modeIndex = getDataEntryLabelIndex(sim[topologyKey][structureKey],"batchId")
            for d in sim[topologyKey][structureKey]["data"]:
                if d[modeIndex] > modeOffset:
                    modeOffset = d[modeIndex]
            modeOffset += 1
        else:
            logger.error(f"Appending mode {mode} not implemented")
            raise Exception("Appending mode not implemented")

    # 1) structure NOT in sim and structure NOT in sim2app
    if not structureInSim and not structureInSim2app:
        #Do nothing
        pass
    # 2) structure NOT in sim and structure in sim2app
    elif not structureInSim and structureInSim2app:
        if topologyKey not in sim:
            sim[topologyKey] = {}
        sim[topologyKey][structureKey] = copy.deepcopy(sim2app[topologyKey][structureKey])

        idOffset   = 0
        modeOffset = 0
    # 3) structure in sim and structure NOT in sim2app
    elif structureInSim and not structureInSim2app:
        #Do nothing
        pass
    # 4) structure in sim and structure in sim2app
    elif structureInSim and structureInSim2app:

        #Handle labels
        simLabels     = sim[topologyKey][structureKey]["labels"]
        sim2appLabels = sim2app[topologyKey][structureKey]["labels"]

        availLabels = ["id","type","resId","chainId","modelId","batchId"]

        for label in simLabels:
            if label not in availLabels:
                logger.error(f"Label {label} not available for structure.")
                raise Exception("Label not available for structure.")

        for label in sim2appLabels:
            if label not in availLabels:
                logger.error(f"Label {label} not available for structure.")
                raise Exception("Label not available for structure.")

        for label in availLabels:
            labelInSim     = (label in simLabels)
            labelInSim2app = (label in sim2appLabels)

            if not labelInSim and not labelInSim2app:
                #Do nothing
                pass
            elif not labelInSim and labelInSim2app:
                sim[topologyKey][structureKey]["labels"].append(label)
                for d in sim[topologyKey][structureKey]["data"]:
                    d.append(0)
            elif labelInSim and not labelInSim2app:
                sim2app[topologyKey][structureKey]["labels"].append(label)
                for d in sim2app[topologyKey][structureKey]["data"]:
                    d.append(0)
            elif labelInSim and labelInSim2app:
                #Do nothing
                pass

        #Data merging
        for structInfo in sim2app[topologyKey][structureKey]["data"]:
            d = []
            for label in sim[topologyKey][structureKey]["labels"]:
                if label == "id":
                    d.append(structInfo[sim2appLabels.index(label)] + idOffset)
                elif label == mode:
                    d.append(structInfo[sim2appLabels.index(label)] + modeOffset)
                else:
                    d.append(structInfo[sim2appLabels.index(label)])
            sim[topologyKey][structureKey]["data"].append(d)

    #Structure appended
    #########################

    return idOffset,modeOffset

