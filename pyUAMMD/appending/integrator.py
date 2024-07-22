import copy

import logging

from deepdiff import DeepDiff

from ..utils.update.groups import updateComponentSetGroupsLists
from ..utils.update.groups import appendGroups

from ..utils.update.ids import updateComponentSetIds

def appendIntegrator(topologyKey,
                     structureKey,
                     integratorKey,
                     sim,sim2app,mode,idOffset,modeOffset,availGroupTypes):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending integrator
    #Integrator structure:
    #   integratorKey: {
    #       "groups": {...} #Groups list, optional
    #       "integrator1": {
    #           "type": ["Class1","SubClass1"],
    #           "parameters": {
    #               "param1": "value1",
    #               "param2": "value2",
    #               ...
    #           }
    #       },
    #       "integrator2": {
    #           "type": ["Class2","SubClass2"],
    #           "parameters": {
    #               "param1": "value1",
    #               "param2": "value2",
    #               ...
    #           }
    #       },
    #       ...,
    #       "schedule": {
    #           "type": ["Schedule",integratorKey],
    #           "labels": ["order",integratorKey,"step"],
    #           "data": [[...,...,...],
    #                    [...,...,...],
    #                    ...]
    #       }
    #   }

    # Four cases:
    # 1) integrator NOT in sim and integrator NOT in sim2app
    # 2) integrator NOT in sim and integrator in sim2app
    # 3) integrator in sim and integrator NOT in sim2app
    # 4) integrator in sim and integrator in sim2app

    integratorInSim     = (integratorKey in sim)
    integratorInSim2app = (integratorKey in sim2app)

    structureInSim     = (topologyKey in sim)
    if structureInSim:
        structureInSim     = (structureKey in sim[topologyKey])

    structureInSim2app = (topologyKey in sim2app)
    if structureInSim2app:
        structureInSim2app    = (structureKey in sim2app[topologyKey])

    if (integratorInSim or integratorInSim2app) and not structureInSim:
        logger.warning("Structure must be present in at least one simulation")
        logger.error("Cannot merge integrator without any structure")
        raise Exception("Cannot merge integrator without structure")

    # 1) integrator NOT in sim and integrator NOT in sim2app
    if not integratorInSim and not integratorInSim2app:
        #Do nothing
        pass
    # 2) integrator NOT in sim and integrator in sim2app
    elif not integratorInSim and integratorInSim2app:
        #In this case the only thing to do is update the sim2app[integratorKey] groups list entries (if any)
        sim[integratorKey] = updateComponentSetGroupsLists(mode,sim2app[integratorKey],idOffset,modeOffset,availGroupTypes)
    # 3) integrator in sim and integrator NOT in sim2app
    elif integratorInSim and not integratorInSim2app:
        #Do nothing
        pass
    # 4) integrator in sim and integrator in sim2app
    elif integratorInSim and integratorInSim2app:
        #We append (after updating) the groupsLists in sim2app[integratorKey] to sim[integratorKey]
        appendGroups(mode,sim[integratorKey],sim2app[integratorKey],idOffset,modeOffset,availGroupTypes)

        integratorsSim     = [k for k in sim[integratorKey].keys() if "Schedule" not in sim[integratorKey][k]["type"]]
        integratorsSim2app = [k for k in sim2app[integratorKey].keys() if "Schedule" not in sim2app[integratorKey][k]["type"]]

        for inte2app in integratorsSim2app:
            if inte2app not in integratorsSim:
                sim[integratorKey][inte2app] = copy.deepcopy(sim2app[integratorKey][inte2app])
            else:
                #If are equal do nothing. Else raise error
                equal = not DeepDiff(sim[integratorKey][inte2app],
                                     sim2app[integratorKey][inte2app],
                                     ignore_order=True,report_repetition=True)
                if not equal:
                    logger.error("Cannot merge integrator with different parameters")
                    raise Exception("Cannot merge integrator with different parameters")

        #Merge schedule
        scheduleName = [k for k in sim[integratorKey].keys() if "Schedule" in sim[integratorKey][k]["type"]]
        if len(scheduleName) == 0:
            logger.error("Cannot merge integrator without schedule")
            raise Exception("Cannot merge integrator without schedule")
        elif len(scheduleName) > 1:
            logger.error("Cannot merge integrator with more than one schedule")
            raise Exception("Cannot merge integrator with more than one schedule")

        scheduleName = scheduleName[0]

        scheduleOffset = 0
        for inte in sim[integratorKey][scheduleName]["data"]:
            if inte[0] > scheduleOffset: #it assumes that the schedule is ordered and well formed
                scheduleOffset = inte[0]

        scheduleName2app = [k for k in sim2app[integratorKey].keys() if "Schedule" in sim2app[integratorKey][k]["type"]]
        if len(scheduleName2app) == 0:
            logger.error("Cannot merge integrator without schedule")
            raise Exception("Cannot merge integrator without schedule")
        elif len(scheduleName2app) > 1:
            logger.error("Cannot merge integrator with more than one schedule")
            raise Exception("Cannot merge integrator with more than one schedule")

        scheduleName2app = scheduleName2app[0]

        for inte in sim2app[integratorKey][scheduleName2app]["data"]:
            inteName = inte[1]
            if inteName not in integratorsSim:
                inte[0] += scheduleOffset
                sim[integratorKey][scheduleName]["data"].append(inte)

    #Integrator appended
    #########################
