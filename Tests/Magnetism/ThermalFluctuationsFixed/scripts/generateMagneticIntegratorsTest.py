import sys,os

import pyUAMMD

import math

import json
import jsbeautifier

with open("parameters.json", "r") as f:
    param = json.load(f)

N = param["N"]

timeStep = param["timeStep"]

temperature = param["temperature"]
volume = param["volume"]

anisotropy = param["anisotropy"]
gyroRatio = param["gyroRatio"]
damping = param["damping"]
msat = param["msat"]
magneticMoment = volume*msat

nSteps        = param["nSteps"]
firstStepSave = param["firstStepSave"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

MIA = param["magneticAlgorithm"] #Magnetic integration algorithm
#Compute box size
L = param["L"]
box = [L,L,L]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "MagneticThermalFluctuationsTest"

simulation["global"] = {}

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]  = [["A", 2, 0.0, 1]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]


simulation["integrator"] = {}
simulation["integrator"]["eulerMaruyamaRigid"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["type"] = ["Magnetic", "Fixed"]
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["msat"] = msat
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["damping"] = damping
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["gyroRatio"] = gyroRatio
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["magneticIntegrationAlgorithm"] = MIA

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyamaRigid", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction", "anisotropy", "magnetization"]
simulation["state"]["data"] = []
for i in range(N):
    simulation["state"]["data"].append([i, [0,0,0], [1.0 ,0 ,0, 0 ], anisotropy, [0,0,1,magneticMoment]])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []
for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

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
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"] = "magnet"
simulation["simulationStep"]["write"]["parameters"]["startStep"] = firstStepSave
"""
simulation["simulationStep"]["eulerMaruyamaRigid"] = {}
simulation["simulationStep"]["eulerMaruyamaRigid"]["parameters"] = {}

"""
#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
