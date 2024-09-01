import os
import numpy as np

import pyUAMMD

nBeads    = 30
nPolymers = 100

sigma = 1.0
L = nBeads*sigma+10.0*sigma

timeStep = 0.001
frictionConstant = 1.0

nSteps = 100000*20000

nStepsInfo   = 1000
nStepsOutput = 10000

Kb = 100.0
Ka = 100.0

compressionVelocity = -0.01

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
simulation["topology"]["structure"]["labels"] = ["id", "type", "modelId"]
simulation["topology"]["structure"]["data"] = []

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads)):
        simulation["state"]["data"].append([particleId, [0.0, 0.0, sigma*j-nBeads*sigma/2.0]])
        simulation["topology"]["structure"]["data"].append([particleId, "A", i])
        particleId += 1

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["bonds"] = {}

simulation["topology"]["forceField"]["bonds"]["type"]       = ["Bond2", "Harmonic"]
simulation["topology"]["forceField"]["bonds"]["parameters"] = {}
simulation["topology"]["forceField"]["bonds"]["labels"]     = ["id_i", "id_j", "K", "r0"]
simulation["topology"]["forceField"]["bonds"]["data"]       = []

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads)-1):
        simulation["topology"]["forceField"]["bonds"]["data"].append([particleId, particleId+1, Kb, sigma])
        particleId += 1
    particleId += 1

simulation["topology"]["forceField"]["angles"] = {}

simulation["topology"]["forceField"]["angles"]["type"]       = ["Bond3", "KratkyPorod"]
simulation["topology"]["forceField"]["angles"]["parameters"] = {}
simulation["topology"]["forceField"]["angles"]["labels"]     = ["id_i", "id_j", "id_k", "K"]
simulation["topology"]["forceField"]["angles"]["data"]       = []

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads)-2):
        simulation["topology"]["forceField"]["angles"]["data"].append([particleId, particleId+1, particleId+2, Ka])
        particleId += 1
    particleId += 2


simulation["topology"]["forceField"]["shell"] = {}
simulation["topology"]["forceField"]["shell"]["type"] = ["External","SphericalShell"]
simulation["topology"]["forceField"]["shell"]["parameters"] = {}
simulation["topology"]["forceField"]["shell"]["parameters"]["shellCenter"] = [0.0, 0.0, 0.0]
simulation["topology"]["forceField"]["shell"]["parameters"]["shellRadius"] = nBeads*sigma+2.0*sigma
simulation["topology"]["forceField"]["shell"]["parameters"]["shellEpsilon"] = 1.0
simulation["topology"]["forceField"]["shell"]["parameters"]["shellSigma"]   = sigma
simulation["topology"]["forceField"]["shell"]["parameters"]["minShellRadius"] = nBeads*sigma/2.0/2.0
simulation["topology"]["forceField"]["shell"]["parameters"]["radiusVelocity"] = compressionVelocity

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
