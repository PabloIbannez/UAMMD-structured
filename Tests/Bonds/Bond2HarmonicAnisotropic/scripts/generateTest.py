import sys,os

import pyUAMMD

import math
import numpy as np
import json
import jsbeautifier

#Import the json input file
with open("parameters.json", "r") as f:
    param = json.load(f)

#Read the parameters
N              = param["N"]

pos0 = 0.0
pos1 = param["r0"]

L    = param["L"]
k    = param["k"]
r0   = param["r0"]
box         = [L,L,L]

temperature = param["temperature"]
nSteps = param["nSteps"]
nStepsOutput = param["nStepsOutput"]
firstStepSave =	param["firstStepSave"]

#Create simulation

simulation = pyUAMMD.simulation()
simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "test"

simulation["global"] = {}
simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "radius", "mass", "charge"]
simulation["global"]["types"]["data"]  = [["A", 1, 0, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}
simulation["integrator"]["eulerMaruyama"] = {}
simulation["integrator"]["eulerMaruyama"]["type"] = ["Brownian", "EulerMaruyama"]
simulation["integrator"]["eulerMaruyama"]["parameters"] = {}
simulation["integrator"]["eulerMaruyama"]["parameters"]["timeStep"] = 0.0001
simulation["integrator"]["eulerMaruyama"]["parameters"]["viscosity"] = 1

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyama", nSteps], #En steps tenfo que poner el número de pasos total, no sé cómo
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []

for i in range(N):
    if i % 2 == 0:
        simulation["state"]["data"].append([i, [pos0, 0, 0]])
    if i % 2 == 1:
        simulation["state"]["data"].append([i, [pos1, 0, 0]])
        
simulation["topology"] = {}
simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["Bond"] = {}
simulation["topology"]["forceField"]["Bond"]["labels"] = ["id_i", "id_j", "K", "r0"]
simulation["topology"]["forceField"]["Bond"]["data"] = []

K = [1*k, 0.5*k, 2.5*k]
R0 = [r0, r0, r0]

for i in range(N):
    if i % 2 ==	0:
        simulation["topology"]["forceField"]["Bond"]["data"].append([i, i + 1, K, R0])

simulation["topology"]["forceField"]["Bond"]["parameters"] = {}
simulation["topology"]["forceField"]["Bond"]["type"] = ["Bond2", "HarmonicAnisotropic"]


simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []

for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

# #Output

simulation["simulationStep"] = {}

simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"]   = "sp"
simulation["simulationStep"]["write"]["parameters"]["startStep"]      = 0

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")
simulation.write("./results/test.json")
