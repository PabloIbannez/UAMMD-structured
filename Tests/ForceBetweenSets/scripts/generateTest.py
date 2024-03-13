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

temperature    = 1

pos0 = param["pos0"]
pos1 = param["pos1"]
L    = param["L"]
k    = param["k"]
r0   = param["r0"]
nPairs = param["nPairs"]

box         = [L,L,L]

intervalStep = 1
nsteps       = 1
dt           = 1e-3
#Create simulation

simulation = pyUAMMD.simulation()
simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "testHessian"

simulation["global"] = {}
simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "radius", "mass", "charge"]
simulation["global"]["types"]["data"]  = [["A", 1, 0, 0], ["B", 1, 0, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}
simulation["integrator"]["BD"] = {}
simulation["integrator"]["BD"]["type"] = ["Brownian", "EulerMaruyama"]
simulation["integrator"]["BD"]["parameters"] = {}
simulation["integrator"]["BD"]["parameters"]["timeStep"] = dt
simulation["integrator"]["BD"]["parameters"]["viscosity"] = 1

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "BD", nsteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []
for i in range(nPairs):
    simulation["state"]["data"].append([2*i,   pos0])
    simulation["state"]["data"].append([2*i+1, pos1])
                                        
simulation["topology"] = {}
simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["Bond"] = {}
simulation["topology"]["forceField"]["Bond"]["labels"] = ["id_i", "id_j", "K", "r0"]
simulation["topology"]["forceField"]["Bond"]["data"] = []
for i in range(nPairs):
    simulation["topology"]["forceField"]["Bond"]["data"].append([2*i, 2*i+1, k, r0])

simulation["topology"]["forceField"]["Bond"]["parameters"] = {}
simulation["topology"]["forceField"]["Bond"]["type"] = ["Bond2", "Harmonic"]

simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []

for i in range(nPairs):
    simulation["topology"]["structure"]["data"].append([2*i,    "A"])
    simulation["topology"]["structure"]["data"].append([2*i+1 , "B"])


# #Output

AIds = [2*i     for i in range(nPairs)]
BIds = [2*i + 1 for i in range(nPairs)]

simulation["simulationStep"] = {}

name = "setsforce"
simulation["simulationStep"][name] = {}
simulation["simulationStep"][name]["type"] = ["MechanicalMeasure", "ForceBetweenSetsMeasure"]
simulation["simulationStep"][name]["parameters"] = {}
simulation["simulationStep"][name]["parameters"]["outputFilePath"]     = "SetsForces.out"
simulation["simulationStep"][name]["parameters"]["intervalStep"]       = intervalStep
simulation["simulationStep"][name]["data"]               = [["1", AIds],["2", BIds]]
simulation["simulationStep"][name]["labels"]             = ["name", "id_list"]

name = "pairforces"
simulation["simulationStep"][name] = {}
simulation["simulationStep"][name]["type"] = ["MechanicalMeasure", "PairwiseForceMeasure"]
simulation["simulationStep"][name]["parameters"] = {}
simulation["simulationStep"][name]["parameters"]["outputFilePath"]     = "PairForces.out"
simulation["simulationStep"][name]["parameters"]["intervalStep"]       = intervalStep
simulation["simulationStep"][name]["parameters"]["mode"]               = "Pairwise_force"


simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"]   = intervalStep
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "pos"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"]   = "sp"

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")
simulation.write("./results/test.json")
