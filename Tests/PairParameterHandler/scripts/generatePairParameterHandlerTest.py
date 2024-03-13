import sys,os

import pyUAMMD

import numpy as np

import json

#######################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

units = param["units"]

N  = param["N"]
N_eps_sgm = param["N_eps_sgm"]

epsilon_min = param["epsilon_min"]
epsilon_max = param["epsilon_max"]
sigma_min   = param["sigma_min"]
sigma_max   = param["sigma_max"]

b  = param["b"]

timeStep    = param["timeStep"]
temperature = param["temperature"]

massA   = param["massA"]
radiusA = param["radiusA"]
chargeA = param["chargeA"]

massB   = param["massB"]
radiusB = param["radiusB"]
chargeB = param["chargeB"]

massC   = param["massC"]
radiusC = param["radiusC"]
chargeC = param["chargeC"]

frictionConstant = param["frictionConstant"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

#Compute box size
L = 10.0*b
box = [L,L,L]

boundParams = param["bond2bound"]
tpy = boundParams["type"]

nonBonded = param["nonBonded"]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "pairParamHandlerTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units",units]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]  = [["A", massA, radiusA, chargeA],
                                          ["B", massB, radiusB, chargeB],
                                          ["C", massC, radiusC, chargeC]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}

simulation["integrator"]["bbk"] = {}
simulation["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["bbk"]["parameters"] = {}
simulation["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["bbk"]["parameters"]["frictionConstant"] = frictionConstant

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "bbk", nSteps]
]

integrators = ["bbk"]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type", "batchId"]
simulation["topology"]["structure"]["data"] = []

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["nl"] = {}
simulation["topology"]["forceField"]["nl"]["type"] = ["VerletConditionalListSet","all"]
simulation["topology"]["forceField"]["nl"]["parameters"] = {"cutOffVerletFactor": 1.2}

simulation["topology"]["forceField"]["Bounds"] = {}
simulation["topology"]["forceField"]["Bounds"]["type"] = ["Bond2", tpy]
simulation["topology"]["forceField"]["Bounds"]["parameters"] = {}
simulation["topology"]["forceField"]["Bounds"]["labels"] = ["id_i", "id_j","maxDistance","K"]
simulation["topology"]["forceField"]["Bounds"]["data"] = []

simulation["topology"]["forceField"]["NonBonded"] = {}
simulation["topology"]["forceField"]["NonBonded"]["type"] = ["NonBonded", nonBonded]
simulation["topology"]["forceField"]["NonBonded"]["parameters"] = {"cutOffFactor":2.5,"condition":"all"}
simulation["topology"]["forceField"]["NonBonded"]["labels"] = ["name_i","name_j","epsilon","sigma","batchId"]
simulation["topology"]["forceField"]["NonBonded"]["data"] = []

simulation["simulationStep"] = {}

partId = 0
simId  = 0
for i in range(N_eps_sgm):

    #Random epsilon sigma
    epsilon = np.random.uniform(epsilon_min, epsilon_max)
    sigma   = np.random.uniform(sigma_min, sigma_max)

    pairs = []
    for n in range(N):
        tpy1 = np.random.choice(["A","B","C"])
        simulation["state"]["data"].append([partId,[0.0,0.0,0.0]])
        simulation["topology"]["structure"]["data"].append([partId, tpy1, simId])
        pairs.append([partId,partId+1])
        partId+=1
        tpy2 = np.random.choice(["A","B","C"])
        simulation["state"]["data"].append([partId,[b,0.0,0.0]])
        simulation["topology"]["structure"]["data"].append([partId, tpy2, simId])
        partId+=1

        simulation["topology"]["forceField"]["NonBonded"]["data"].append([tpy1,tpy2,epsilon,sigma,simId])

        simId+=1

    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"] = {}
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["type"] = ["ParticlesListMeasure", "DistancesMeasure"]
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["parameters"] = {}
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["parameters"]["intervalStep"] = nStepsOutput
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["parameters"]["outputFilePath"] = f"distances_{epsilon}_{sigma}.dat"
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["labels"] = ["id_i", "id_j"]
    simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["data"] = []

    for p in pairs:
        simulation["simulationStep"][f"distance_{epsilon}_{sigma}"]["data"].append(p)

    boundParams = param["bond2bound"]

    K        = boundParams["K"]
    distance = boundParams["distance"]

    for p in pairs:
        simulation["topology"]["forceField"]["Bounds"]["data"].append(p + [distance,K])

#Output

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")


