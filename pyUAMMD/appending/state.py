import copy

import logging

from ..utils.common import getDataEntryLabelIndex

def appendState(topologyKey,
                structureKey,
                stateKey,
                sim,sim2app,mode,idOffset,modeOffset):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending state
    #State structure:
    #   stateKey: {
    #       "labels": ["id","position",...]
    #       "data": [[...],
    #                [...],
    #                ...]

    # Four cases:
    # 1) state NOT in sim and state NOT in sim2app
    # 2) state NOT in sim and state in sim2app
    # 3) state in sim and state NOT in sim2app
    # 4) state in sim and state in sim2app

    stateInSim     = (stateKey in sim)
    stateInSim2app = (stateKey in sim2app)

    structureInSim     = (topologyKey in sim)
    if structureInSim:
        structureInSim     = (structureKey in sim[topologyKey])

    if (stateInSim or stateInSim2app) and not structureInSim:
        logger.warning("Structure must be present in at least one simulation")
        logger.error("Cannot merge state without any structure")
        raise Exception("Cannot merge state without structure")

    # 1) state NOT in sim and state NOT in sim2app
    if not stateInSim and not stateInSim2app:
        #Do nothing
        pass
    # 2) state NOT in sim and state in sim2app
    elif not stateInSim and stateInSim2app:
        sim[stateKey] = copy.deepcopy(sim2app[stateKey])
    # 3) state in sim and state NOT in sim2app
    elif stateInSim and not stateInSim2app:
        #Do nothing
        pass
    # 4) state in sim and state in sim2app
    elif stateInSim and stateInSim2app:
        #At this point, idOffset should be defined from the structure merging

        #Handle labels
        simLabels     = sim[stateKey]["labels"]
        sim2appLabels = sim2app[stateKey]["labels"]

        availLabels = ["position","velocity","direction"]
        availLabelZeros = {"position":[0.0,0.0,0.0],
                           "velocity":[0.0,0.0,0.0],
                           "direction":[0.0,0.0,0.0,1.0]}

        if "id" not in simLabels or "id" not in sim2appLabels:
            logger.error("Cannot merge state without id")
            raise Exception("Cannot merge state without id")

        for label in availLabels:
            labelInSim     = (label in simLabels)
            labelInSim2app = (label in sim2appLabels)

            if not labelInSim and not labelInSim2app:
                #Do nothing
                pass
            elif not labelInSim and labelInSim2app:
                sim[stateKey]["labels"].append(label)
                for d in sim[stateKey]["data"]:
                    d.append(availLabelZeros[label])
            elif labelInSim and not labelInSim2app:
                sim2app[stateKey]["labels"].append(label)
                for d in sim2app[stateKey]["data"]:
                    d.append(availLabelZeros[label])
            elif labelInSim and labelInSim2app:
                #Do nothing
                pass

        #Data merging
        for stateInfo in sim2app[stateKey]["data"]:
            d = []
            for label in sim[stateKey]["labels"]:
                if label == "id":
                    d.append(stateInfo[sim2appLabels.index(label)] + idOffset)
                else:
                    d.append(stateInfo[sim2appLabels.index(label)])
            sim[stateKey]["data"].append(d)

    #State appended
    #########################

