import copy

import logging

from deepdiff import DeepDiff

def appendSystem(systemKey,sim,sim2app,mode):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending system
    #System structure:
    #   systemKey: {
    #       "information": {
    #           "type": ["Simulation","Information"],
    #           "parameters": {"name": simName,
    #                          ...}
    #       }
    #       "backup": {
    #           "type": ["Simulation","Backup"],
    #           "parameters": {...},
    #       }

    # Four cases:
    # 1) system NOT in sim and system NOT in sim2app
    # 2) system NOT in sim and system in sim2app
    # 3) system in sim and system NOT in sim2app
    # 4) system in sim and system in sim2app

    systemInSim     = (systemKey in sim)
    systemInSim2app = (systemKey in sim2app)

    # 1) system NOT in sim and system NOT in sim2app
    if not systemInSim and not systemInSim2app:
        #Do nothing
        pass
    # 2) system NOT in sim and system in sim2app
    elif not systemInSim and systemInSim2app:
        sim[systemKey] = copy.deepcopy(sim2app[systemKey])
    # 3) system in sim and system NOT in sim2app
    elif systemInSim and not systemInSim2app:
        #Do nothing
        pass
    # 4) system in sim and system in sim2app
    elif systemInSim and systemInSim2app:
        infoInSim   = False
        backupInSim = False
        for key in sim[systemKey]:
            tpy,subTpy = sim[systemKey][key]["type"]

            if tpy == "Simulation" and subTpy == "Information":
                if not infoInSim:
                    infoInSim    = True
                    infoInSimKey = key
                else:
                    logger.error("More than one information system in simulation!")
                    raise Exception("More than one information system in simulation!")
            elif tpy == "Simulation" and subTpy == "Backup":
                if not backupInSim:
                    backupInSim    = True
                    backupInSimKey = key
                else:
                    logger.error("More than one backup system in simulation!")
                    raise Exception("More than one backup system in simulation!")
            else:
                logger.error("Unknown system type in simulation!")
                raise Exception("Unknown system type!")

        infoInSim2app   = False
        backupInSim2app = False
        for key in sim2app[systemKey]:
            tpy,subTpy = sim2app[systemKey][key]["type"]

            if tpy == "Simulation" and subTpy == "Information":
                if not infoInSim2app:
                    infoInSim2app    = True
                    infoInSim2appKey = key
                else:
                    logger.error("More than one information system in simulation!")
                    raise Exception("More than one information system in simulation!")
            elif tpy == "Simulation" and subTpy == "Backup":
                if not backupInSim2app:
                    backupInSim2app    = True
                    backupInSim2appKey = key
                else:
                    logger.error("More than one backup system in simulation!")
                    raise Exception("More than one backup system in simulation!")
            else:
                logger.error("Unknown system type in simulation to append!")
                raise Exception("Unknown system type!")

        # Four cases for info:
        # 1) info NOT in sim and info NOT in sim2app
        # 2) info NOT in sim and info in sim2app
        # 3) info in sim and info NOT in sim2app
        # 4) info in sim and info in sim2app

        # 1) info NOT in sim and info NOT in sim2app
        if not infoInSim and not infoInSim2app:
            #Do nothing
            pass
        # 2) info NOT in sim and info in sim2app
        elif not infoInSim and infoInSim2app:
            sim[systemKey][infoInSim2appKey] = copy.deepcopy(sim2app[systemKey][infoInSim2appKey])
        # 3) info in sim and info NOT in sim2app
        elif infoInSim and not infoInSim2app:
            #Do nothing
            pass
        # 4) info in sim and info in sim2app
        elif infoInSim and infoInSim2app:

            nameInSim   = sim[systemKey][infoInSimKey]["parameters"].get("name",None)
            nameInSim2app = sim2app[systemKey][infoInSim2appKey]["parameters"].get("name",None)

            if nameInSim == None and nameInSim2app == None:
                #Do nothing
                pass
            elif nameInSim != None and nameInSim2app == None:
                #Do nothing
                pass
            elif nameInSim == None and nameInSim2app != None:
                sim[systemKey][infoInSimKey]["parameters"]["name"] = nameInSim2app
            elif nameInSim != None and nameInSim2app != None:
                if mode == "modelId":
                    #Name must be the same
                    self.logger.error("System name must be the same. (Mode: modelId)")
                    raise Exception("System name must be the same.")
                elif mode == "batchId":
                    #Combine names
                    sim[systemKey][infoInSimKey]["parameters"]["name"] += "_" + nameInSim2app

            ###################################################################################

            seedInSim     = sim[systemKey][infoInSimKey]["parameters"].get("seed",None)
            seedInSim2app = sim2app[systemKey][infoInSim2appKey]["parameters"].get("seed",None)

            if seedInSim == None and seedInSim2app == None:
                #Do nothing
                pass
            elif seedInSim != None and seedInSim2app == None:
                #Do nothing
                pass
            elif seedInSim == None and seedInSim2app != None:
                sim[systemKey][infoInSimKey]["parameters"]["seed"] = seedInSim2app
            elif seedInSim != None and seedInSim2app != None:
                if seedInSim != seedInSim2app:
                    logger.error("Seeds of two systems are different!")
                    raise Exception("Seeds of two systems are different!")

                sim[systemKey][infoInSimKey]["parameters"]["seed"] = seedInSim2app

        # Four cases for backup:
        # 1) backup NOT in sim and backup NOT in sim2app
        # 2) backup NOT in sim and backup in sim2app
        # 3) backup in sim and backup NOT in sim2app
        # 4) backup in sim and backup in sim2app

        # 1) backup NOT in sim and backup NOT in sim2app
        if not backupInSim and not backupInSim2app:
            #Do nothing
            pass
        # 2) backup NOT in sim and backup in sim2app
        elif not backupInSim and backupInSim2app:
            sim[systemKey][backupInSim2appKey] = copy.deepcopy(sim2app[systemKey][backupInSim2appKey])
        # 3) backup in sim and backup NOT in sim2app
        elif backupInSim and not backupInSim2app:
            #Do nothing
            pass
        # 4) backup in sim and backup in sim2app
        elif backupInSim and backupInSim2app:

            equal = not DeepDiff(sim[systemKey][backupInSimKey],
                                 sim2app[systemKey][backupInSim2appKey],
                                 ignore_order=True,report_repetition=True)

            if not equal:
                logger.error("Backup systems are different!")
                raise Exception("Backup systems are different!")

    #System appended
    #########################

