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

radius    = param["radius"]
viscosity = param["viscosity"]
tolerance = param["tolerance"]
vwall     = param["vwall"]
omega     = param["omega"]
z0        = param["z0"]
H         = param["H"]
k         = param["k"]
Lxy       = param["Lxy"]
mass_p    = param["particleMass"]
mass_f    = param["fluidMass"]
overtone  = param["overtone"]
memory    = param["memory"]
damping   = param["damping"]

toleranceConvergence = param["toleranceConvergence"]
notAcceleratedInterval = param["notAcceleratedInterval"]


f0 = omega/(2*np.pi*overtone)
particleDensity = mass_p/(2*radius)**3
fluidDensity    = mass_f/(2*radius)**3

iterationsList = param["iterationsList"]
nSteps        = 1
timeStep      = 0
nStepsOutput  = 1
nStepsMeasure = 1
temperature   = 0

box = [Lxy,Lxy,2*H]
#Create simulation

simulation = pyUAMMD.simulation()


simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "VQCM"

simulation["global"] = {}
simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "radius", "mass", "charge"]
simulation["global"]["types"]["data"]  = [["A", radius, mass_p, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}
simulation["integrator"]["BBK"] = {}
simulation["integrator"]["BBK"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["BBK"]["parameters"] = {}
simulation["integrator"]["BBK"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["BBK"]["parameters"]["frictionConstant"] = 1

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "BBK", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = [[0, [0,0,-H+z0]]]


simulation["topology"] = {}

simulation["topology"] = {}

simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["Bond"] = {}
simulation["topology"]["forceField"]["Bond"]["labels"] = ["id_i", "K", "r0", "position"]
simulation["topology"]["forceField"]["Bond"]["data"] = [[0, [k, k, k], [0,0,0], [0,0,-H+z0]]]

simulation["topology"]["forceField"]["Bond"]["parameters"] = {}
simulation["topology"]["forceField"]["Bond"]["type"] = ["Bond1", "FixedHarmonic"]

simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = [[0, "A"]]



#Output

simulation["simulationStep"] = {}

for nIterations in iterationsList:
    name = "vqcmb"+str(nIterations)
    
    simulation["simulationStep"][name] = {}
    simulation["simulationStep"][name]["type"] = ["VQCMMeasure", "VQCMMeasure"]
    
    simulation["simulationStep"][name]["parameters"] = {}
    simulation["simulationStep"][name]["parameters"]["kernel"]          = "Gaussian"
    simulation["simulationStep"][name]["parameters"]["intervalStep"]         = nStepsOutput
    simulation["simulationStep"][name]["parameters"]["tolerance"]            = tolerance
    simulation["simulationStep"][name]["parameters"]["f0"]                   = f0
    simulation["simulationStep"][name]["parameters"]["overtone"]             = [overtone]
    simulation["simulationStep"][name]["parameters"]["hydrodynamicRadius"]   = radius
    simulation["simulationStep"][name]["parameters"]["viscosity"]            = viscosity
    simulation["simulationStep"][name]["parameters"]["vwall"]                = vwall
    simulation["simulationStep"][name]["parameters"]["fluidDensity"]         = fluidDensity
    simulation["simulationStep"][name]["parameters"]["maxNIterations"]       = nIterations
    simulation["simulationStep"][name]["parameters"]["toleranceConvergence"] = toleranceConvergence
    simulation["simulationStep"][name]["parameters"]["maxNIterations"]         = nIterations
    simulation["simulationStep"][name]["parameters"]["memory"]                 = memory
    simulation["simulationStep"][name]["parameters"]["damping"]                = damping
    simulation["simulationStep"][name]["parameters"]["notAcceleratedInterval"] = notAcceleratedInterval
    simulation["simulationStep"][name]["labels"] = ["ids", "positions", "step"]
    simulation["simulationStep"][name]["data"] = [[0, [0,0,-H+z0],0]]


    #simulation["simulationStep"][name]["parameters"]["nTethers"]           = nTethers
    
    filename = "Impedanceb"+str(nIterations)+".out"
    simulation["simulationStep"][name]["parameters"]["outputFilePath"]     = filename


simulation["simulationStep"]["vqcmMob"] = {}
simulation["simulationStep"]["vqcmMob"]["type"] = ["VQCMMeasure", "VQCMMeasureFromMobility"]

simulation["simulationStep"]["vqcmMob"]["parameters"] = {}
simulation["simulationStep"]["vqcmMob"]["parameters"]["intervalStep"]       = nStepsOutput
simulation["simulationStep"]["vqcmMob"]["parameters"]["tolerance"]          = tolerance
simulation["simulationStep"]["vqcmMob"]["parameters"]["omega"]              = omega
simulation["simulationStep"]["vqcmMob"]["parameters"]["hydrodynamicRadius"] = radius
simulation["simulationStep"]["vqcmMob"]["parameters"]["viscosity"]          = viscosity
simulation["simulationStep"]["vqcmMob"]["parameters"]["vwall"]              = vwall
simulation["simulationStep"]["vqcmMob"]["parameters"]["fluidDensity"]       = fluidDensity


filename = "ImpedanceFromMobility.out"
simulation["simulationStep"]["vqcmMob"]["parameters"]["outputFilePath"]     = filename

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")
simulation.write("./results/simulationOneParticle.json")
