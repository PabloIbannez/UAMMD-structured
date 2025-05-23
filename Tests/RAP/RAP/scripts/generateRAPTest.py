import sys,os

import pyUAMMD

import numpy as np

import json
import jsbeautifier

#######################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

N  = param["N"]
b  = param["b"]

timeStep = param["timeStep"]

temperature = param["temperature"]

mass   = param["mass"]
radius = param["radius"]

Krap = param["Krap"]
R    = param["R"]

viscosity = param["viscosity"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]

#Compute box size
L = 2.0*b*(N-1)
box = [L,L,L]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "RAPTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","None"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}

simulation["integrator"]["eulerMaruyamaRigid"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["type"] = ["Brownian", "EulerMaruyamaRigidBody"]
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["viscosity"] = viscosity

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyamaRigid", nSteps]
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction"]
simulation["state"]["data"] = []

previous = np.zeros(3)
for i in range(N):
    # Generate a random quaternion
    q = np.random.rand(4)
    q /= np.linalg.norm(q)
    q = q.tolist()
    simulation["state"]["data"].append([i,[0.0,0.0,0.0],q])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []

for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["rap"] = {}
simulation["topology"]["forceField"]["rap"]["type"]   = ["Bond2", "RAP"]
simulation["topology"]["forceField"]["rap"]["labels"] = ["id_i", "id_j", "K", "R"]
simulation["topology"]["forceField"]["rap"]["data"]   = []

for i in range(N-1):
    simulation["topology"]["forceField"]["rap"]["data"].append([i, i+1, Krap, R])

#Output
simulation["simulationStep"] = {}
simulation["simulationStep"][f"write"] = {}
simulation["simulationStep"][f"write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"][f"write"]["parameters"] = {}
simulation["simulationStep"][f"write"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"][f"write"]["parameters"]["outputFilePath"] = f"output"
simulation["simulationStep"][f"write"]["parameters"]["outputFormat"]   = "itpd"

simulation["simulationStep"][f"writeSPO"] = {}
simulation["simulationStep"][f"writeSPO"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"][f"writeSPO"]["parameters"] = {}
simulation["simulationStep"][f"writeSPO"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"][f"writeSPO"]["parameters"]["outputFilePath"] = f"output"
simulation["simulationStep"][f"writeSPO"]["parameters"]["outputFormat"]   = "spo"

simulation["simulationStep"]["thermo"] = {}
simulation["simulationStep"]["thermo"]["type"] = ["ThermodynamicMeasure","ThermodynamicQuantityMeasure"]
simulation["simulationStep"]["thermo"]["parameters"] = {}
simulation["simulationStep"]["thermo"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["thermo"]["parameters"]["outputFilePath"] = "thermo.dat"

#simulation["simulationStep"]["potMeasure"] = {}
#simulation["simulationStep"]["potMeasure"]["type"] = ["ParticlesListMeasure","PotentialMeasure"]
#simulation["simulationStep"]["potMeasure"]["parameters"] = {
#        "intervalStep": 1,
#        "outputFilePath": "pot.dat"
#        }
#simulation["simulationStep"]["potMeasure"]["labels"] = ["id"]
#simulation["simulationStep"]["potMeasure"]["data"] = [[i] for i in range(N)]

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
