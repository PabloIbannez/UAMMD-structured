import sys,os

import pyUAMMD

import numpy as np

import json
import jsbeautifier

#######################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

units = "None"

N  = 8000
#Compute box size
box = [25,25,25]

timeStep         = param["timeStep"]
frictionConstant = param["frictionConstant"]

temperature = 1.0

epsilon  = 1.0
sigma    = 1.0
rcFactor = 2.5

mass   = 1.0
radius = np.power(2.0,1.0/6.0)*sigma/2.0

nSteps        = param["nSteps"]
nStepsInfo    = param["nStepsInfo"]
nStepsOutput  = param["nStepsOutput"]

nModel = param["nModel"]
nSim   = param["nSim"]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "conditionalListSetTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units",units]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}

simulation["integrator"]["bbk"] = {}
simulation["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["bbk"]["parameters"] = {}
simulation["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["bbk"]["parameters"]["frictionConstant"] = frictionConstant

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "bbk", nSteps]
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type", "modelId", "batchId"]
simulation["topology"]["structure"]["data"] = []

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["groups"] = {}
simulation["topology"]["forceField"]["groups"]["type"] = ["Groups", "GroupsList"]
simulation["topology"]["forceField"]["groups"]["parameters"] = {}
simulation["topology"]["forceField"]["groups"]["labels"] = ["name","type","selection"]
simulation["topology"]["forceField"]["groups"]["data"] = []

simulation["topology"]["forceField"]["verletList"] = {}
simulation["topology"]["forceField"]["verletList"]["type"] = ["VerletConditionalListSet", "intra_inter"]
simulation["topology"]["forceField"]["verletList"]["parameters"] = {}
simulation["topology"]["forceField"]["verletList"]["parameters"]["cutOffVerletFactor"] = 1.5

simulation["simulationStep"] = {}

simulation["simulationStep"] = {}
simulation["simulationStep"]["groups"] = {}
simulation["simulationStep"]["groups"]["type"] = ["Groups", "GroupsList"]
simulation["simulationStep"]["groups"]["parameters"] = {}
simulation["simulationStep"]["groups"]["labels"] = ["name","type","selection"]
simulation["simulationStep"]["groups"]["data"] = []

index  = 0
offset = 0;
for s in range(nSim):
    for mdl in range(nModel):

        #Create particles
        n       = int(np.ceil(np.power(N,1.0/3.0)))
        spacingX = box[0]/n
        spacingY = box[1]/n
        spacingZ = box[2]/n

        print("[INFO] Creating particles ..., spacing =", spacingX, spacingY, spacingZ,"diameter =", 2.0*radius)

        positions = []
        while len(positions) < N:
            for x in range(n):
                for y in range(n):
                    for z in range(n):
                        positions.append([x*spacingX,y*spacingY,z*spacingZ])

        for i,p in enumerate(positions):
            simulation["state"]["data"].append([i+offset, p])
            simulation["topology"]["structure"]["data"].append([i+offset, "A", mdl, s])

        simulation["topology"]["forceField"]["groups"]["data"].append(["pg"+str(index), "ModelIdsBatchIds", [[mdl,s]]])

        simulation["topology"]["forceField"][f"lj{index}"] = {}
        simulation["topology"]["forceField"][f"lj{index}"]["type"] = ["NonBonded", "LennardJonesType1"]
        simulation["topology"]["forceField"][f"lj{index}"]["parameters"] = {"cutOffFactor": rcFactor,
                                                                            "condition":"intra",
                                                                            "group":"pg"+str(index)}
        simulation["topology"]["forceField"][f"lj{index}"]["labels"] = ["name_i", "name_j", "epsilon", "sigma"]
        simulation["topology"]["forceField"][f"lj{index}"]["data"] = [["A", "A", epsilon, sigma]]

        simulation["simulationStep"]["groups"]["data"].append(["pg"+str(index), "ModelIdsBatchIds", [[mdl,s]]])

        simulation["simulationStep"][f"write{index}"] = {}
        simulation["simulationStep"][f"write{index}"]["type"] = ["WriteStep", "WriteStep"]
        simulation["simulationStep"][f"write{index}"]["parameters"] = {"group" : f"pg{index}",
                                                                       "startStep" : nStepsOutput}
        simulation["simulationStep"][f"write{index}"]["parameters"]["intervalStep"] = nStepsOutput
        simulation["simulationStep"][f"write{index}"]["parameters"]["outputFilePath"] = f"output{index}"
        simulation["simulationStep"][f"write{index}"]["parameters"]["outputFormat"] = "sp"


        index  += 1
        offset += N

#Output

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsInfo

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")


