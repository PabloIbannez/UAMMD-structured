import os
import numpy as np

import pyUAMMD

NpartPerBatch  = 2e4
Nbatch         = 100

conc = 0.1

L = (NpartPerBatch / conc)**(1.0/3.0)

sigma = 1.0

timeStep = 0.001
frictionConstant = 1.0

nSteps = 100000*20000

nStepsInfo   = 1000
nStepsOutput = 10000

#########################

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "LennardJones"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","None"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", 1.0, sigma/2.0, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[[L,L,L], 1.0]]

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
simulation["topology"]["structure"]["labels"] = ["id", "type", "batchId"]
simulation["topology"]["structure"]["data"] = []

Linterval = np.linspace(-(L-2.0*sigma)/2.0, (L-2.0*sigma)/2.0, int(np.ceil(np.cbrt(NpartPerBatch))+1))

# Convert Linterval to cyclic array
Linterval = np.roll(Linterval, int(np.ceil(np.cbrt(NpartPerBatch))/2.0))

particleId = 0
while particleId < NpartPerBatch*Nbatch:
    for lx in Linterval:
        for ly in Linterval:
            for lz in Linterval:

                batchId = int(particleId/NpartPerBatch)

                if particleId < NpartPerBatch*Nbatch:
                    simulation["state"]["data"].append([particleId, [lx, ly, lz]])
                    simulation["topology"]["structure"]["data"].append([particleId, "A", batchId])
                    particleId += 1



simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["nl"] = {}
simulation["topology"]["forceField"]["nl"]["type"] = ["VerletConditionalListSet","all"]
simulation["topology"]["forceField"]["nl"]["parameters"] = {"cutOffVerletFactor": 1.2}

simulation["topology"]["forceField"]["lj"] = {}
simulation["topology"]["forceField"]["lj"]["type"] = ["NonBonded", "LennardJonesType2"]
simulation["topology"]["forceField"]["lj"]["parameters"] = {}
simulation["topology"]["forceField"]["lj"]["parameters"]["condition"]    = "all"
simulation["topology"]["forceField"]["lj"]["parameters"]["cutOffFactor"] = 2.5

simulation["topology"]["forceField"]["lj"]["labels"] = ["name_i","name_j","epsilon","sigma"]
simulation["topology"]["forceField"]["lj"]["data"] = [["A","A",1.0,sigma]]

simulation["simulationStep"] = {}

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsInfo

simulation["simulationStep"]["output"] = {}
simulation["simulationStep"]["output"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["output"]["parameters"] = {"intervalStep": nStepsOutput,
                                                        "outputFilePath": "output",
                                                        "outputFormat": "sp"}

#Check if ./simulation folder exists, if not create it
if not os.path.exists("./simulation"):
    os.makedirs("./simulation")

print("Writing simulation file...")
simulation.write("./simulation/simulation.json")
