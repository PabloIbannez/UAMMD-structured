import os,sys
import uuid

from typing import final

import copy
import itertools

import json

import logging

from deepdiff import DeepDiff

from .utils.common import getDataEntryLabelIndex

from .utils.update.ids import updateDataEntryIds
from .utils.update.ids import updateComponentSetIds

from .utils.update.groups import getGroupsListEntriesFromComponentSet
from .utils.update.groups import updateComponentSetGroupsLists
from .utils.update.groups import appendGroups

from .appending import appendSystem
from .appending import appendGlobal
from .appending import appendStructure
from .appending import appendState
from .appending import appendIntegrator
from .appending import appendSimulationStep
from .appending import appendForceField
from .appending import appendPatchyParticles

class simulation:

    availModes      = ["batchId","modelId"]

    availGroupTypes = ["Ids","notIds",
                       "Types","notTypes",
                       "ModelIds","notModelIds",
                       "BatchIds","notBatchIds",
                       "ModelIdsBatchIds","notModelIdsBatchIds"]

    id_labels = ["id",
                 "id_i","id_j","id_k","id_l"]

    id_list_labels = ["id_list",
                      "idSet_i","idSet_j","idSet_k","idSet_l"]

    type_labels = ["name",
                   "name_i","name_j"]

    def __init__(self,simulationJSON=None,debug=False):

        self.localID = uuid.uuid1().hex

        ########################

        self.logger = logging.getLogger("pyUAMMD")

        if debug:
            self.logger.setLevel(logging.DEBUG) #<----
        else:
            self.logger.setLevel(logging.INFO) #<----

        self.logger.debug("Initializing simulation, ID: {}".format(self.localID))

        ########################

        if simulationJSON is not None:
            self.sim = copy.deepcopy(simulationJSON)
        else:
            self.sim = {}

        ########################

        #Checks

        stateInSim     = ("state" in self.sim)
        topologyInSim  = ("topology" in self.sim)

        if topologyInSim:
            structureInSim = ("structure" in self.sim["topology"])
        else:
            structureInSim = False

        if stateInSim and not structureInSim:
            self.logger.error("Added state but not structure")
            raise Exception("Added state but not structure")

        if stateInSim and structureInSim:
            Nsim = len(self.sim["state"]["data"])
            Nstr = len(self.sim["topology"]["structure"]["data"])

            if Nsim != Nstr:
                self.logger.error("Number of particles in state and structure does not match")
                raise Exception("Number of particles in state and structure does not match")

        ########################

        self.logger.debug("Simulation created")


    def getID(self):
        return self.localID

    def getNumberOfParticles(self):
        try:
            return len(self.sim["topology"]["structure"]["data"])
        except:
            return 0

    def append(self,sim2app,mode: str = "batchId") -> None:

        sim2appID = sim2app.getID()

        # Check if ID are the same
        if self.localID == sim2appID:
            self.logger.error("Cannot append simulation to itself")
            raise Exception("Cannot append simulation to itself")

        self.logger.debug(f"Appending \"{sim2appID}\" to \"{self.localID}\" with mode \"{mode}\"")
        if mode not in self.availModes:
            self.logger.error(f"Trying append with mode \"{mode}\", but this mode is not available. " \
                             f"Available modes are: {self.availModes}")
            raise Exception(f"Trying append with not available mode \"{mode}\"")

        if mode == "batchId":
            #Check both simulations have exactly the same entries.
            #In principle this is not compulsory but avoid a large
            #number of errors.

            def get_keys(dictionary):
                result = []
                for key, value in dictionary.items():
                    if type(value) is dict:
                        new_keys = get_keys(value)
                        result.append(key)
                        for innerkey in new_keys:
                            result.append(f'{key}/{innerkey}')
                    else:
                        result.append(key)
                return result

            simEntries     = get_keys(self.sim)
            sim2appEntries = get_keys(sim2app.sim)

            #Remove all simulationStep entries

            simEntries = [entry for entry in simEntries if "simulationStep" not in entry]
            sim2appEntries = [entry for entry in sim2appEntries if "simulationStep" not in entry]

            if len(simEntries) != len(sim2appEntries):
                self.logger.error("Trying to append simulations with different entries while using mode \"batchId\"")
                raise Exception("Trying to append simulations with different entries")
            else:
                if not (sorted(simEntries) == sorted(sim2appEntries)):
                    #Compute not shared entries
                    self.logger.error("Trying to append simulations with different entries (but same length) while using mode \"batchId\"")
                    raise Exception("Trying to append simulations with different entries")

        ########################################

        #For avoiding problems patchy particles are appended first
        #We need reference sim structure to be unmodified
        appendPatchyParticles("topology",
                              "structure",
                              "forceField",
                              self,sim2app,mode,
                              self.availGroupTypes,
                              self.id_labels, self.id_list_labels,self.type_labels)

        ###########################################################


        appendSystem("system",self,sim2app,mode)
        appendGlobal("global",self,sim2app,mode)

        idOffset,modeOffset = appendStructure("topology","structure",self,sim2app,mode)

        appendState("topology",
                    "structure",
                    "state",self,sim2app,mode,idOffset,modeOffset)
        appendIntegrator("topology",
                         "structure",
                         "integrator",self,sim2app,mode,idOffset,modeOffset,self.availGroupTypes)

        appendSimulationStep("topology",
                             "structure",
                             "simulationStep",
                             self,sim2app,mode,idOffset,modeOffset,self.availGroupTypes,
                             self.id_labels, self.id_list_labels)
        appendForceField("topology",
                         "structure",
                         "forceField",
                         self,sim2app,mode,idOffset,modeOffset,self.availGroupTypes,
                         self.id_labels, self.id_list_labels,self.type_labels,
                         ignoredEntriesType=["PatchyParticles"])

    def write(self,output:str,legacy:bool=False):
        self.logger.debug("Writing to json...")

        try:
            if legacy:
                raise Exception("Forcing legacy mode")
            import JFIO
            self.logger.debug("Writing with JFIO")
            JFIO.write(output,self.sim,formatted=True)
        except:
            import jsbeautifier
            self.logger.warning("Writing with legacy")
            with open(output,"w") as f:
                opts = jsbeautifier.default_options()
                opts.indent_size  = 2
                f.write(jsbeautifier.beautify(json.dumps(self.sim), opts))

    def setValue(self,path,value):
        simTmp = self.sim
        for p in path[:-1]:
            simTmp = simTmp[p]
        simTmp[path[-1]] = value

    #__eq__ is used in many places, so it is better to define it
    #Do not touch it !!!
    @final
    def __eq__(self, other):
        return id(self) == id(other)

    #With the following methods simulation behaves like a dict

    def __getitem__(self,key):
        return self.sim[key]

    def __setitem__(self,key,value):
        self.sim[key]=value

    def __delitem__(self, key):
        del self.sim[key]

    def __contains__(self, key):
        return key in self.sim

    def keys(self):
        return list(self.sim.keys())

    def values(self):
        return list(self.sim.values())

    def items(self):
        return list(self.sim.items())

    def pop(self, key, default=None):
        try:
            value = self.sim[key]
            del self.sim[key]
            return value
        except KeyError:
            if default is not None:
                return default
            else:
                raise

    def run(self):
        self.logger.debug("Runing simulation...")

        try:
            from .utils.launcher.pyUAMMDlauncher import UAMMDlauncher
            UAMMDlauncher(self.sim)
        except:
            self.logger.error("Something went wrong with UAMMDlauncher")
