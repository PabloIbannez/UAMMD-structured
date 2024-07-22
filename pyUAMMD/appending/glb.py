import copy
import logging

from deepdiff import DeepDiff

from ..utils.common import getDataEntryLabelIndex

def appendGlobal(globalKey,sim,sim2app,mode):

    logger = logging.getLogger("pyUAMMD")

    #########################
    #Appending global data
    #Global data structure:
    #   globalKey: {
    #       "units": {
    #           "type": ["Units","UnitsName"],
    #       },
    #       "types": {
    #           "type": ["Types","TypesName"],
    #       },
    #       "fundamentals": {
    #           "type": ["Fundamentals","FundamentalsName"],
    #           "parameters": {"param1": value1,
    #                         "param2": value2, ...},
    #       },
    #       "ensemble": {
    #           "type": ["Ensemble","EnsembleName"],
    #           "labels": ["ensembleParam1","ensembleParam2",...],
    #           "data": [[value1],[value2],...]
    #       }
    #   }

    # Four cases:
    # 1) global NOT in sim and global NOT in sim2app
    # 2) global NOT in sim and global in sim2app
    # 3) global in sim and global NOT in sim2app
    # 4) global in sim and global in sim2app

    globalInSim     = (globalKey in sim)
    globalInSim2app = (globalKey in sim2app)

    # 1) global NOT in sim and global NOT in sim2app
    if not globalInSim and not globalInSim2app:
        #Do nothing
        pass
    # 2) global NOT in sim and global in sim2app
    elif not globalInSim and globalInSim2app:
        sim[globalKey] = copy.deepcopy(sim2app[globalKey])
    # 3) global in sim and global NOT in sim2app
    elif globalInSim and not globalInSim2app:
        #Do nothing
        pass
    # 4) global in sim and global in sim2app
    elif globalInSim and globalInSim2app:

        #Process parameters

        unitsInSim     = ("units" in sim[globalKey])
        unitsInSim2app = ("units" in sim2app[globalKey])

        typesInSim     = ("types" in sim[globalKey])
        typesInSim2app = ("types" in sim2app[globalKey])

        fundamentalsInSim     = ("fundamentals" in sim[globalKey])
        fundamentalsInSim2app = ("fundamentals" in sim2app[globalKey])

        ensembleInSim     = ("ensemble" in sim[globalKey])
        ensembleInSim2app = ("ensemble" in sim2app[globalKey])

        #Units

        # 1) units NOT in sim and units NOT in sim2app
        if not unitsInSim and not unitsInSim2app:
            #Do nothing
            pass
        # 2) units NOT in sim and units in sim2app
        elif not unitsInSim and unitsInSim2app:
            sim[globalKey]["units"] = copy.deepcopy(sim2app[globalKey]["units"])
        # 3) units in sim and units NOT in sim2app
        elif unitsInSim and not unitsInSim2app:
            #Do nothing
            pass
        # 4) units in sim and units in sim2app
        elif unitsInSim and unitsInSim2app:
            #Check if units are the same
            equal = not DeepDiff(sim[globalKey]["units"],sim2app[globalKey]["units"])
            if not equal:
                logger.error("Units must match between simulation and appending")
                raise Exception("Units do not match")
            else:
                #Do nothing
                pass

        #Types

        # 1) types NOT in sim and types NOT in sim2app
        if not typesInSim and not typesInSim2app:
            #Do nothing
            pass
        # 2) types NOT in sim and types in sim2app
        elif not typesInSim and typesInSim2app:
            sim[globalKey]["types"] = copy.deepcopy(sim2app[globalKey]["types"])
        # 3) types in sim and types NOT in sim2app
        elif typesInSim and not typesInSim2app:
            #Do nothing
            pass
        # 4) types in sim and types in sim2app
        elif typesInSim and typesInSim2app:
            #Check type is the same
            if sim[globalKey]["types"]["type"] != sim2app[globalKey]["types"]["type"]:
                logger.error("Type for types must match between simulation and appending")
                raise Exception("Types do not match")

            #Check label is the same
            if sim[globalKey]["types"]["labels"] != sim2app[globalKey]["types"]["labels"]:
                logger.error("Labels for types must match between simulation and appending")
                raise Exception("Types labels do not match")

            #Iterate over data in sim2app, if data is not in sim, add it
            for d in sim2app[globalKey]["types"]["data"]:
                if d not in sim[globalKey]["types"]["data"]:
                    sim[globalKey]["types"]["data"].append(d)

        #Fundamentals

        # 1) fundamentals NOT in sim and fundamentals NOT in sim2app
        if not fundamentalsInSim and not fundamentalsInSim2app:
            #Do nothing
            pass
        # 2) fundamentals NOT in sim and fundamentals in sim2app
        elif not fundamentalsInSim and fundamentalsInSim2app:
            sim[globalKey]["fundamentals"] = copy.deepcopy(sim2app[globalKey]["fundamentals"])
        # 3) fundamentals in sim and fundamentals NOT in sim2app
        elif fundamentalsInSim and not fundamentalsInSim2app:
            #Do nothing
            pass
        # 4) fundamentals in sim and fundamentals in sim2app
        elif fundamentalsInSim and fundamentalsInSim2app:
            #Check type is the same
            if sim[globalKey]["fundamentals"]["type"] != sim2app[globalKey]["fundamentals"]["type"]:
                logger.error("Type for fundamentals must match between simulation and appending")
                raise Exception("Fundamentals do not match")

            #Check parameters are the same
            if sim[globalKey]["fundamentals"]["parameters"] != sim2app[globalKey]["fundamentals"]["parameters"]:
                logger.error("Parameters for fundamentals must match between simulation and appending")
                raise Exception("Fundamentals parameters do not match")

        #Ensemble

        # 1) ensemble NOT in sim and ensemble NOT in sim2app
        if not ensembleInSim and not ensembleInSim2app:
            #Do nothing
            pass
        # 2) ensemble NOT in sim and ensemble in sim2app
        elif not ensembleInSim and ensembleInSim2app:
            sim[globalKey]["ensemble"] = copy.deepcopy(sim2app[globalKey]["ensemble"])
        # 3) ensemble in sim and ensemble NOT in sim2app
        elif ensembleInSim and not ensembleInSim2app:
            #Do nothing
            pass
        # 4) ensemble in sim and ensemble in sim2app
        elif ensembleInSim and ensembleInSim2app:
            #Check type is the same
            if sim[globalKey]["ensemble"]["type"] != sim2app[globalKey]["ensemble"]["type"]:
                logger.error("Type for ensemble must match between simulation and appending")
                raise Exception("Ensemble do not match")

            #Check label is the same
            if sim[globalKey]["ensemble"]["labels"] != sim2app[globalKey]["ensemble"]["labels"]:
                logger.error("Labels for ensemble must match between simulation and appending")
                raise Exception("Ensemble labels do not match")

            #Check data is the same
            if sim[globalKey]["ensemble"]["data"] != sim2app[globalKey]["ensemble"]["data"]:
                logger.error("Data for ensemble must match between simulation and appending")
                raise Exception("Ensemble data do not match")

    #Global data appended
    #########################
