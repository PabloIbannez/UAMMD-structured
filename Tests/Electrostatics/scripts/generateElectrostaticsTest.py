import os
import numpy as np

import pyUAMMD

N    = 2e4
conc = 0.1

L = (N / conc)**(1.0/3.0)

sigma = 1.0
chg   = 1.0

timeStep = 0.02
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
simulation["system"]["info"]["parameters"]["name"] = "Electrostatics"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","KcalMol_A"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", 1.0, sigma/2.0, chg]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[[L,L,L], 1.0]]

simulation["integrator"] = {}

simulation["integrator"]["bbk"] = {}
simulation["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["bbk"]["parameters"] = {}
simulation["integrator"]["bbk"]["parameters"]["timeStep"]         = timeStep
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

Linterval = np.linspace(-(L-2.0*sigma)/2.0, (L-2.0*sigma)/2.0, int(np.ceil(np.cbrt(N))+1))

particleId = 0
while particleId < N:
    for lx in Linterval:
        for ly in Linterval:
            for lz in Linterval:

                if particleId < N:
                    simulation["state"]["data"].append([particleId, [lx, ly, lz]])
                    particleId += 1

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []

for i in range(int(N)):
    simulation["topology"]["structure"]["data"].append([i, "A"])

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["elect"] = {}
simulation["topology"]["forceField"]["elect"]["type"] = ["LongRange", "Electrostatics"]
simulation["topology"]["forceField"]["elect"]["parameters"] = {}
simulation["topology"]["forceField"]["elect"]["parameters"]["gaussianWidth"]      = sigma
simulation["topology"]["forceField"]["elect"]["parameters"]["dielectricConstant"] = 1.0
simulation["topology"]["forceField"]["elect"]["parameters"]["tolerance"]          = 1e-4

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
if not os.path.exists("./results"):
    os.makedirs("./results")

print("Writing simulation file ...")
simulation.write("./results/simulation.json")
