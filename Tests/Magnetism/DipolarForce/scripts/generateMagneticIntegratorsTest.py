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
N              = 1
temperature    = 0
permeability   = param["permeability"]
volume         = param["volume"]
msat           = param["msat"]
viscosity      = param["viscosity"]
rij            = param["rij"]
magneticMoment = msat*volume
radius         = (3*volume/(4*math.pi))**(1./3.)

nSteps        = param["nSteps"]
timeStep      = param["timeStep"]
firstStepSave = param["firstStepSave"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

MIA = "none"
#Compute box size


dist = (rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])**0.5
L = 20*dist
box = [L,L,L]
#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "Dipolar_interactions"

simulation["global"] = {}

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "radius", "mass", "charge"]
simulation["global"]["types"]["data"]  = [["A", radius, 0, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}
simulation["integrator"]["magneticFixed"] = {}
simulation["integrator"]["magneticFixed"]["type"] = ["Magnetic", "Brownian"]
simulation["integrator"]["magneticFixed"]["parameters"] = {}
simulation["integrator"]["magneticFixed"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["magneticFixed"]["parameters"]["msat"] = msat
simulation["integrator"]["magneticFixed"]["parameters"]["viscosity"] = viscosity
simulation["integrator"]["magneticFixed"]["parameters"]["magneticIntegrationAlgorithm"] = MIA

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "magneticFixed", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction", "magnetization"]
simulation["state"]["data"] = []

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type","modelId"]
simulation["topology"]["structure"]["data"] = []

simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["nl"] = {}
simulation["topology"]["forceField"]["nl"]["type"] = ["VerletConditionalListSet","intra_inter"]
simulation["topology"]["forceField"]["nl"]["parameters"] = {"cutOffVerletFactor": 1.1}


simulation["topology"]["forceField"]["NonBonded"] = {}
simulation["topology"]["forceField"]["NonBonded"]["type"]       = ["NonBonded", "DipolarMagnetic"]
simulation["topology"]["forceField"]["NonBonded"]["parameters"] = {"cutOff":1.5*dist,
                                                                   "condition":"intra",
                                                                   "permeability": permeability}


simulation["state"]["data"].append([0, [0,0,0], [1.0 ,0 ,0, 0 ], [0,0,1,magneticMoment],0.0])
simulation["state"]["data"].append([1, rij, [1.0 , 0, 0, 0], [0, 0, 1, magneticMoment],0.0])

simulation["topology"]["structure"]["data"].append([0,   "A", 0])
simulation["topology"]["structure"]["data"].append([1, "A", 0])

#Output

simulation["simulationStep"] = {}
simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"]   = "spm"
simulation["simulationStep"]["write"]["parameters"]["startStep"]      = firstStepSave

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
