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
anisotropy  = param["anisotropy"]
gyroRatio   = param["gyroRatio"]
damping     = param["damping"]
viscosity   = param["viscosity"]


magneticMoment = 4.*math.pi/3.*msat*radius**3

nSteps        = param["nSteps"]
timeStep      = param["timeStep"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

#Compute box size
L = param["L"]
box = [L,L,L]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "RelaxometryNeel"

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
simulation["integrator"]["LLG-Brown"] = {}
simulation["integrator"]["LLG-Brown"]["type"] = ["Magnetic", "Brownian"]
simulation["integrator"]["LLG-Brown"]["parameters"] = {}
simulation["integrator"]["LLG-Brown"]["parameters"]["timeStep"]  = timeStep
simulation["integrator"]["LLG-Brown"]["parameters"]["msat"]      = msat
simulation["integrator"]["LLG-Brown"]["parameters"]["damping"]   = damping
simulation["integrator"]["LLG-Brown"]["parameters"]["viscosity"] = viscosity
simulation["integrator"]["LLG-Brown"]["parameters"]["gyroRatio"] = gyroRatio
simulation["integrator"]["LLG-Brown"]["parameters"]["magneticIntegrationAlgorithm"] = "LLG_Heun"

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "LLG-Brown", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction", "magnetization", "anisotropy"]
simulation["state"]["data"] = []
for i in range(N):
    simulation["state"]["data"].append([i, [0,0,0], [1.0 ,0 ,0, 0 ],
                                        [0,0,1,magneticMoment], anisotropy])

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
simulation["simulationStep"]["write"]["type"] = ["MagneticMeasure", "MeasureMeanMagnetization"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "RelaxationCombined.magnet"

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulationCombined.json")
