import logging

import copy
import itertools

from tqdm import tqdm

def mergeSimulationsSet(simulationsSet,
                        independentSimulationSteps=True,
                        additionalIgnoredSimulationStepTypes=[]):

    ignoredSimulationStepTypes = ['UtilsStep','FlowControl','Groups'] + additionalIgnoredSimulationStepTypes

    logger = logging.getLogger("pyUAMMD")

    if independentSimulationSteps:
        #We assume that batchId is 0 (not set) for each independent simulation
        #Check it
        for sim in simulationsSet:
            structureLabels = sim["topology"]["structure"]["labels"]
            if "batchId" in structureLabels:
                logger.error("[pyUAMMD] BatchId already set, for a single simulation. Cannot aggregate simulations")
                raise Exception("BatchId already set")

        #Update simulationsStep for each simulation in the simulation set
        #This ensures the simulation step is applied over particles for the
        #same batchId
        for simIndex,sim in enumerate(simulationsSet):

            if "simulationStep" in sim.keys():
                groupDefinitionRequired = False

                keysToRename = []
                for simStep in sim["simulationStep"].keys():

                    simStepType = sim["simulationStep"][simStep]["type"][0]

                    isIgnored = False
                    for ignoredType in ignoredSimulationStepTypes:
                        if ignoredType in simStepType:
                            isIgnored = True
                            break

                    if not isIgnored:

                        #Each simulationStep apply to a batchId group
                        if "group" not in sim["simulationStep"][simStep]["parameters"].keys():
                            sim["simulationStep"][simStep]["parameters"]["group"] = f"batchId_{simIndex}"
                            groupDefinitionRequired = True

                        keysToRename.append(simStep)

                    if "Groups" in simStepType:

                        #We rename all groups declared
                        for i in range(len(sim["simulationStep"][simStep]["data"])):
                            oldName = sim["simulationStep"][simStep]["data"][i][0]
                            newName = oldName + f"_{simIndex}"
                            sim["simulationStep"][simStep]["data"][i][0] = newName

                            #Rename all references to the group
                            for simStep2rename in sim["simulationStep"].keys():
                                if "group" in sim["simulationStep"][simStep2rename]["parameters"].keys():
                                    if sim["simulationStep"][simStep2rename]["parameters"]["group"] == oldName:
                                        sim["simulationStep"][simStep2rename]["parameters"]["group"] = newName

                for k in keysToRename:
                    sim["simulationStep"][k+f"_{simIndex}"] = sim["simulationStep"].pop(k)

                if groupDefinitionRequired:
                    if "groups_batchId" not in sim["simulationStep"].keys():
                        sim["simulationStep"]["groups_batchId"] = {
                            "type": ["Groups","GroupsList"],
                            "parameters": {},
                            "labels":["name","type","selection"],
                            "data":[
                                [f"batchId_{simIndex}","BatchIds",[0]],
                            ]
                        }
                    else:
                        logger.error("[pyUAMMD] groups_batchId already defined")
                        raise Exception("groups_batchId already defined")

    aggregatedSimulation = None
    for sim in tqdm(simulationsSet,desc="Merging simulations set"):
        if aggregatedSimulation is None:
            aggregatedSimulation = copy.deepcopy(sim)
        else:
            aggregatedSimulation.append(sim,mode="batchId")

    #Check if more than one neighbor list is defined
    #List all neighbor lists
    neighborLists = []
    if "forceField" in aggregatedSimulation["topology"].keys():
        for entry in aggregatedSimulation["topology"]["forceField"].keys():
            if aggregatedSimulation["topology"]["forceField"][entry]["type"][0] == "VerletConditionalListSet":
                neighborLists.append(entry)

    if len(neighborLists) > 1:
        nlsTypes = ""
        for i in neighborLists:
            nlsTypes += f"\"{i}\" ({aggregatedSimulation['topology']['forceField'][i]['type'][1]}), "
        nlsTypes = nlsTypes[:-2]+"."
        logger.warning(f"[pyUAMMD] Resulsting simulation will have multiple neighbor lists !!!. Neighbor lists types: {nlsTypes}")

    return aggregatedSimulation

