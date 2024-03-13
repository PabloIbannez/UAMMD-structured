import sys,os

import pyUAMMD

import numpy as np

import json

from tqdm import tqdm

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

b = param["b"]

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

boundParams = param["bound"]
tpy = boundParams["type"]
boundEps = boundParams["epsilon"]
boundSig = boundParams["sigma"]
boundSurfacePos = boundParams["surfacePosition"]

surface = param["surface"]

#Create simulation

sim = pyUAMMD.simulation()

sim["system"] = {}
sim["system"]["info"] = {}
sim["system"]["info"]["type"] = ["Simulation","Information"]
sim["system"]["info"]["parameters"] = {}
sim["system"]["info"]["parameters"]["name"] = "singleParamHandlerTest"

sim["global"] = {}

sim["global"]["units"] = {}
sim["global"]["units"]["type"] = ["Units",units]

sim["global"]["types"] = {}
sim["global"]["types"]["type"]   = ["Types","Basic"]
sim["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
sim["global"]["types"]["data"]  = [["A", massA, radiusA, chargeA],
                                   ["B", massB, radiusB, chargeB],
                                   ["C", massC, radiusC, chargeC]]

sim["global"]["ensemble"] = {}
sim["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
sim["global"]["ensemble"]["labels"] = ["box", "temperature"]
sim["global"]["ensemble"]["data"]   = [[box, temperature]]

sim["integrator"] = {}

sim["integrator"]["bbk"] = {}
sim["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
sim["integrator"]["bbk"]["parameters"] = {}
sim["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
sim["integrator"]["bbk"]["parameters"]["frictionConstant"] = frictionConstant

sim["integrator"]["schedule"] = {}
sim["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
sim["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
sim["integrator"]["schedule"]["data"] = [
    [1, "bbk", nSteps]
]

integrators = ["bbk"]

sim["state"] = {}
sim["state"]["labels"] = ["id", "position"]
sim["state"]["data"] = []

sim["topology"] = {}
sim["topology"]["structure"] = {}
sim["topology"]["structure"]["labels"] = ["id", "type", "batchId"]
sim["topology"]["structure"]["data"] = []

sim["topology"]["forceField"] = {}

sim["topology"]["forceField"]["Bounds"] = {}
sim["topology"]["forceField"]["Bounds"]["type"] = ["Surface", tpy]
sim["topology"]["forceField"]["Bounds"]["parameters"] = {"surfacePosition": boundSurfacePos}
sim["topology"]["forceField"]["Bounds"]["labels"] = ["name","epsilon","sigma"]
sim["topology"]["forceField"]["Bounds"]["data"] = [
    ["A", boundEps, boundSig],
    ["B", boundEps, boundSig],
    ["C", boundEps, boundSig]
]

sim["topology"]["forceField"]["Surface"] = {}
sim["topology"]["forceField"]["Surface"]["type"] = ["Surface", surface]
sim["topology"]["forceField"]["Surface"]["parameters"] = {}
sim["topology"]["forceField"]["Surface"]["labels"] = ["name","epsilon","sigma","batchId"]
sim["topology"]["forceField"]["Surface"]["data"] = []

sim["simulationStep"] = {}
sim["simulationStep"]["groups"] = {}
sim["simulationStep"]["groups"]["type"] = ["Groups", "GroupsList"]
sim["simulationStep"]["groups"]["parameters"] = {}
sim["simulationStep"]["groups"]["labels"] = ["name","type","selection"]
sim["simulationStep"]["groups"]["data"] = []

for i in range(N_eps_sgm):
    sim["simulationStep"]["groups"]["data"].append([f"pg{i}","BatchIds",[i]])

partId = 0
for i in tqdm(range(N_eps_sgm)):

    #Random epsilon sigma
    epsilon = np.random.uniform(epsilon_min, epsilon_max)
    sigma   = np.random.uniform(sigma_min, sigma_max)

    for n in range(N):
        tpy1 = np.random.choice(["A","B","C"])
        sim["state"]["data"].append([partId,[0.0,0.0,b]])
        sim["topology"]["structure"]["data"].append([partId, tpy1, i])
        partId+=1

    sim["topology"]["forceField"]["Surface"]["data"].append(["A",epsilon,sigma,i])
    sim["topology"]["forceField"]["Surface"]["data"].append(["B",epsilon,sigma,i])
    sim["topology"]["forceField"]["Surface"]["data"].append(["C",epsilon,sigma,i])

    sim["simulationStep"][f"position_{epsilon}_{sigma}"] = {}
    sim["simulationStep"][f"position_{epsilon}_{sigma}"]["type"] = ["WriteStep", "WriteStep"]
    sim["simulationStep"][f"position_{epsilon}_{sigma}"]["parameters"] = {"group":f"pg{i}"}
    sim["simulationStep"][f"position_{epsilon}_{sigma}"]["parameters"]["intervalStep"] = nStepsOutput
    sim["simulationStep"][f"position_{epsilon}_{sigma}"]["parameters"]["outputFilePath"] = f"position_{epsilon}_{sigma}"
    sim["simulationStep"][f"position_{epsilon}_{sigma}"]["parameters"]["outputFormat"] = "sp"

#Output

sim["simulationStep"]["info"] = {}
sim["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
sim["simulationStep"]["info"]["parameters"] = {}
sim["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

sim.write("./results/simulation.json")
