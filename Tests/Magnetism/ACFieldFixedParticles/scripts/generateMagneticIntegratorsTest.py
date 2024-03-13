import sys,os

import pyUAMMD

import math

import json
import jsbeautifier

with open("parameters.json", "r") as f:
    param = json.load(f)

N = param["N"]


temperature = param["temperature"]
radius      = param["radius"]
msat        = param["msat"]
b0          = param["b0"]
anisotropy  = param["anisotropy"]
gyroRatio   = param["gyroRatio"]
damping     = param["damping"]
frequency   = param["frequency"]


nCycles = param["nCycles"]
stepsPerCycle = param["stepsPerCycle"]
nSteps        = nCycles*stepsPerCycle
timeStep      = 1/(stepsPerCycle*frequency)
firstStepSave = stepsPerCycle*(nCycles-1)
nStepsOutput  = stepsPerCycle/param["nOutputPerCycle"]
nStepsMeasure = stepsPerCycle/param["nMeasurePerCycle"]

magneticMoment = 4.*math.pi/3.*msat*radius**3
MIA = "LLG_Heun"

#Compute box size
L = param["L"]
box = [L,L,L]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "ACMagneticField"

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
simulation["integrator"]["magneticFixed"]["type"] = ["Magnetic", "Fixed"]
simulation["integrator"]["magneticFixed"]["parameters"] = {}
simulation["integrator"]["magneticFixed"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["magneticFixed"]["parameters"]["msat"] = msat
simulation["integrator"]["magneticFixed"]["parameters"]["damping"] = damping
simulation["integrator"]["magneticFixed"]["parameters"]["gyroRatio"] = gyroRatio
simulation["integrator"]["magneticFixed"]["parameters"]["magneticIntegrationAlgorithm"] = MIA

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "magneticFixed", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction", "magnetization", "anisotropy"]
simulation["state"]["data"] = []
for i in range(N):
    simulation["state"]["data"].append([i, [0,0,0], [1.0 ,0 ,0, 0 ], [0,0,1,magneticMoment], anisotropy])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []
for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])


simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["External"] = {}
simulation["topology"]["forceField"]["External"]["type"] = ["External", "ACMagneticField"]
simulation["topology"]["forceField"]["External"]["parameters"] = {}
simulation["topology"]["forceField"]["External"]["parameters"]["b0"] = b0
simulation["topology"]["forceField"]["External"]["parameters"]["frequency"] = frequency
simulation["topology"]["forceField"]["External"]["parameters"]["direction"] = [0,0,1]

#Output

simulation["simulationStep"] = {}
simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput


simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["MagneticMeasure", "MeasureMeanMagnetization"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output.magnet"
#simulation["simulationStep"]["write"]["parameters"]["outputFormat"] = "magnet"
simulation["simulationStep"]["write"]["parameters"]["startStep"] = firstStepSave

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
