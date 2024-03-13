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
N              = 4
temperature    = 0

pos0 = param["pos0"]
pos1 = param["pos1"]
pos2 = param["pos2"]
pos3 = param["pos3"]

phi0 = param["phi0"]
n    = param["n"]
L    = param["L"]
k    = param["k"]

box         = [L,L,L]
temperature = 0

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
simulation["global"]["types"]["data"]  = [["A", 1, 0, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}
simulation["integrator"]["BBK"] = {}
simulation["integrator"]["BBK"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["BBK"]["parameters"] = {}
simulation["integrator"]["BBK"]["parameters"]["timeStep"] = 0.001
simulation["integrator"]["BBK"]["parameters"]["frictionConstant"] = 1

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "BBK", 1],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = [[0, pos0],
                               [1, pos1],
                               [2, pos2],
                               [3, pos3]]


simulation["topology"] = {}
simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["Bond"] = {}
simulation["topology"]["forceField"]["Bond"]["labels"] = ["id_i", "id_j", "id_k", "id_l", "K", "n","phi0"]
simulation["topology"]["forceField"]["Bond"]["data"] = [[0, 1, 2, 3, k, n, phi0]]

simulation["topology"]["forceField"]["Bond"]["parameters"] = {}
simulation["topology"]["forceField"]["Bond"]["type"] = ["Bond4", "Dihedral"]

simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = [[0, "A"],
                                               [1, "A"],
                                               [2, "A"],
                                               [3, "A"]]



# #Output

simulation["simulationStep"] = {}

name = "Hessian_analytical"
simulation["simulationStep"][name] = {}
simulation["simulationStep"][name]["type"] = ["MechanicalMeasure", "HessianMeasure"]
simulation["simulationStep"][name]["parameters"] = {}
simulation["simulationStep"][name]["parameters"]["mode"]               = "Analytical"
simulation["simulationStep"][name]["parameters"]["outputFilePath"]     = "HessianAnalytical.out"
simulation["simulationStep"][name]["parameters"]["intervalStep"]       = 1

name = "Hessian_numerical"
simulation["simulationStep"][name] = {}
simulation["simulationStep"][name]["type"] = ["MechanicalMeasure", "HessianMeasure"]
simulation["simulationStep"][name]["parameters"] = {}
simulation["simulationStep"][name]["parameters"]["mode"]               = "Numerical"
simulation["simulationStep"][name]["parameters"]["outputFilePath"]     = "HessianNumerical.out"
simulation["simulationStep"][name]["parameters"]["intervalStep"]       = 1


#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")
simulation.write("./results/test.json")
