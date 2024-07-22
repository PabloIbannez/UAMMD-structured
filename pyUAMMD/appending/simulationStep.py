import copy

from deepdiff import DeepDiff

import logging

from ..utils.update.groups import appendGroups

from ..utils.update.ids import updateComponentSetIds

def appendSimulationStep(topologyKey,
                         structureKey,
                         simulationStepKey,
                         sim,sim2app,mode,idOffset,modeOffset,availGroupTypes,
                         id_labels,id_list_labels):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending simulation step
    #Simulation steps structure:
    #   simulationStepKey: {
    #       "groups": {...} #Groups list, optional
    #       "step1": {
    #           "type": ["Class1","SubClass1"],
    #           "parameters": {
    #               "intervalStep": "value",
    #               "startStep": "value", #optional
    #               "endStep": "value", #optional
    #               "param2": "value2",
    #               ...
    #           }
    #           "labels": ["label1","label2",...],
    #           "data": [[...],
    #                    [...],
    #                    ...]
    #       },
    #       "step2": {
    #           "type": ["Class2","SubClass2"],
    #           "parameters": {
    #               "intervalStep": "value",
    #               "startStep": "value", #optional
    #               "endStep": "value", #optional
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
    # 1) simulationStep NOT in sim and simulationStep NOT in sim2app
    # 2) simulationStep NOT in sim and simulationStep in sim2app
    # 3) simulationStep in sim and simulationStep NOT in sim2app
    # 4) simulationStep in sim and simulationStep in sim2app

    simulationStepInSim     = (simulationStepKey in sim)
    simulationStepInSim2app = (simulationStepKey in sim2app)

    structureInSim     = (topologyKey in sim)
    if structureInSim:
        structureInSim     = (structureKey in sim[topologyKey])

    structureInSim2app = (topologyKey in sim2app)
    if structureInSim2app:
        structureInSim2app    = (structureKey in sim2app[topologyKey])

    if (simulationStepInSim or simulationStepInSim2app) and not structureInSim:
        logger.warning("Structure must be present in at least one simulation")
        logger.error("Cannot merge simulationStep without any structure")
        raise Exception("Cannot merge simulationStep without structure")

    # 1) simulationStep NOT in sim and simulationStep NOT in sim2app
    if not simulationStepInSim and not simulationStepInSim2app:
        #Do nothing
        pass
    # 2) simulationStep NOT in sim and simulationStep in sim2app
    elif not simulationStepInSim and simulationStepInSim2app:
        #Check if the sim2app have particles defined.
        if structureInSim2app:
            #In this case, we first update the groups in sim2app and then de ids in the remainder entries
            sim[simulationStepKey] = updateComponentSetIds(mode,
                                                               updateComponentSetGroupsLists(mode,
                                                                                             sim2app[simulationStepKey],
                                                                                             idOffset,modeOffset,availGroupTypes),
                                                               idOffset,id_labels,id_list_labels)
        else:
            sim[simulationStepKey] = copy.deepcopy(sim2app[simulationStepKey])

    # 3) simulationStep in sim and simulationStep NOT in sim2app
    elif simulationStepInSim and not simulationStepInSim2app:
        #Do nothing
        pass
    # 4) simulationStep in sim and simulationStep in sim2app
    elif simulationStepInSim and simulationStepInSim2app:

        #Check if the sim2app have particles defined.
        if structureInSim2app:
            #We append (after updating) the groupsLists in sim2app[simulationStepKey] to sim[simulationStepKey]
            updatedSimulationStep2app = updateComponentSetIds(mode,sim2app[simulationStepKey],idOffset,id_labels,id_list_labels)
            appendGroups(mode,sim[simulationStepKey],updatedSimulationStep2app,idOffset,modeOffset,availGroupTypes)
        else:
            updatedSimulationStep2app = copy.deepcopy(sim2app[simulationStepKey])
            appendGroups(mode,sim[simulationStepKey],updatedSimulationStep2app,0,0,availGroupTypes)

        simulationStepsSim     = [k for k in sim[simulationStepKey].keys() if sim[simulationStepKey][k]["type"] != ["Groups","GroupsList"]]
        simulationStepsSim2app = [k for k in sim2app[simulationStepKey].keys() if sim2app[simulationStepKey][k]["type"] != ["Groups","GroupsList"]]

        for simStep2app in simulationStepsSim2app:
            if simStep2app not in simulationStepsSim:
                sim[simulationStepKey][simStep2app] = copy.deepcopy(updatedSimulationStep2app[simStep2app])
            else:
                equalParam = not DeepDiff(sim[simulationStepKey][simStep2app],
                                          updatedSimulationStep2app[simStep2app],
                                          ignore_order=True,report_repetition=True,
                                          exclude_paths=["root['data']"]) #Note data is excluded from comparison
                if equalParam:

                    if "data" not in sim[simulationStepKey][simStep2app] and "data" not in updatedSimulationStep2app[simStep2app]:
                        continue

                    dataSim     = sim[simulationStepKey][simStep2app]["data"]
                    dataSim2app = updatedSimulationStep2app[simStep2app]["data"]

                    if len(dataSim) != len(dataSim2app):
                        equalData = False
                    else:
                        equalData = not DeepDiff(dataSim,dataSim2app,ignore_order=True,report_repetition=True)

                    if equalData:
                        pass
                    else:
                        sim[simulationStepKey][simStep2app]["data"].extend(updatedSimulationStep2app[simStep2app]["data"])

                else:
                    logger.error("Only simulation step entries which differ in data can be appended. Simulation step entry: \"{}\"".format(simStep2app))
                    raise Exception("Only simulation step entries which differ in data can be appended")

    #Simulation step appended
    #########################
