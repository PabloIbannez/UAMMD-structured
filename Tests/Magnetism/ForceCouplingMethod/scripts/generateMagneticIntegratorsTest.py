import sys,os

import pyUAMMD

import math

import json
import jsbeautifier

with open("parameters.json", "r") as f:
    param = json.load(f)

N           = 1
temperature = 0
tolerance   = param["tolerance"]
radius      = param["radius"]
force       = param["force"]
viscosity   = param["viscosity"]
lbox        = param["lbox"]

nSteps        = param["nSteps"]
timeStep      = param["timeStep"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]


MIA = "none" 
#Compute box size
box = [lbox, lbox, lbox]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "FCM"

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
simulation["integrator"]["FCM"] = {}
simulation["integrator"]["FCM"]["type"] = ["Magnetic", "ForceCouplingMethod"]
simulation["integrator"]["FCM"]["parameters"] = {}
simulation["integrator"]["FCM"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["FCM"]["parameters"]["msat"] = 1
simulation["integrator"]["FCM"]["parameters"]["viscosity"] = viscosity
simulation["integrator"]["FCM"]["parameters"]["tolerance"] = tolerance
simulation["integrator"]["FCM"]["parameters"]["magneticIntegrationAlgorithm"] = MIA
simulation["integrator"]["FCM"]["parameters"]["hydrodynamicRadius"] = radius

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "FCM", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction"]
simulation["state"]["data"] = []
for i in range(N):
    simulation["state"]["data"].append([i, [0,0,0], [1.0 ,0 ,0, 0 ]])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []
for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])
    
    
simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["External"] = {}
simulation["topology"]["forceField"]["External"]["type"] = ["External", "ConstantForce"]
simulation["topology"]["forceField"]["External"]["parameters"] = {}
simulation["topology"]["forceField"]["External"]["parameters"]["constantForce"] = [0,0,force]

#Output

simulation["simulationStep"] = {}
simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput


simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "pos"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"] = "sp"

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
