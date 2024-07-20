import sys,os

import pyUAMMD

import numpy as np

import json
import jsbeautifier

#######################

#COMPONENTS_PATH = os.path.join(os.getenv('UAMMD_PATH'),
#                               'USCM/Components.json')

COMPONENTS_PATH = os.path.join("../../",
                               'USCM/Components.json')

with open(COMPONENTS_PATH) as f:
    components = json.load(f)

#######################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

units = param["units"]

N  = param["N"]
b  = param["b"]

timeStep = param["timeStep"]

temperature = param["temperature"]
lambd       = param["lambda"]

mass   = param["mass"]
radius = param["radius"]
charge = param["charge"]

frictionConstant = param["frictionConstant"]
viscosity = param["viscosity"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

#Compute box size
L = 10.0*b
box = [L,L,L]

nonBondedList = param["nonBondedList"]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "nonbondedTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units",units]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, charge]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVTlambda"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature","lambda"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature,lambd]]

simulation["integrator"] = {}

simulation["integrator"]["bbk"] = {}
simulation["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
simulation["integrator"]["bbk"]["parameters"] = {}
simulation["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["bbk"]["parameters"]["frictionConstant"] = frictionConstant

simulation["integrator"]["gjf"] = {}
simulation["integrator"]["gjf"]["type"] = ["Langevin", "GJF"]
simulation["integrator"]["gjf"]["parameters"] = {}
simulation["integrator"]["gjf"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["gjf"]["parameters"]["frictionConstant"] = frictionConstant

simulation["integrator"]["eulerMaruyama"] = {}
simulation["integrator"]["eulerMaruyama"]["type"] = ["Brownian", "EulerMaruyama"]
simulation["integrator"]["eulerMaruyama"]["parameters"] = {}
simulation["integrator"]["eulerMaruyama"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["eulerMaruyama"]["parameters"]["viscosity"] = viscosity

simulation["integrator"]["eulerMaruyamaRigid"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["type"] = ["Brownian", "EulerMaruyamaRigidBody"]
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["viscosity"] = viscosity

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyama", nSteps],
    [2, "eulerMaruyamaRigid", nSteps],
    [3, "bbk", nSteps],
    [4, "gjf", nSteps]
]

integrators = ["eulerMaruyama", "eulerMaruyamaRigid", "bbk", "gjf"]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type", "batchId"]
simulation["topology"]["structure"]["data"] = []

simulation["topology"]["forceField"] = {}

boundParams = param["bond2bound"]
tpy = boundParams["type"]

simulation["topology"]["forceField"]["nl"] = {}
simulation["topology"]["forceField"]["nl"]["type"] = ["VerletConditionalListSet","all"]
simulation["topology"]["forceField"]["nl"]["parameters"] = {"cutOffVerletFactor": 1.2}

simulation["topology"]["forceField"]["Bounds"] = {}
simulation["topology"]["forceField"]["Bounds"]["type"] = ["Bond2", tpy]
simulation["topology"]["forceField"]["Bounds"]["parameters"] = {}
simulation["topology"]["forceField"]["Bounds"]["labels"] = ["id_i", "id_j","maxDistance","K"]
simulation["topology"]["forceField"]["Bounds"]["data"] = []

simulation["topology"]["forceField"]["groups"] = {}
simulation["topology"]["forceField"]["groups"]["type"] = ["Groups", "GroupsList"]
simulation["topology"]["forceField"]["groups"]["parameters"] = {}
simulation["topology"]["forceField"]["groups"]["labels"] = ["name","type","selection"]
simulation["topology"]["forceField"]["groups"]["data"] = []

simulation["simulationStep"] = {}

partId = 0
simId  = 0
for nonBnd in components["Interactor"]["Pair"]:

    cls,subCls,f = nonBnd

    if subCls not in nonBondedList:
        print("[WARNING] NonBonded type {} not found".format(subCls))
        continue

    pairs = []
    for i in range(N):
        simulation["state"]["data"].append([partId,[0.0,0.0,0.0]])
        simulation["topology"]["structure"]["data"].append([partId, "A", simId])
        pairs.append([partId,partId+1])
        partId+=1
        simulation["state"]["data"].append([partId,[b,0.0,0.0]])
        simulation["topology"]["structure"]["data"].append([partId, "A", simId])
        partId+=1
        simId+=1

    simulation["topology"]["forceField"]["groups"]["data"].append([f"pg{subCls}","Ids",[x for y in pairs for x in y]])

    simulation["topology"]["forceField"][subCls] = {}
    simulation["topology"]["forceField"][subCls]["type"] = [cls, subCls]
    simulation["topology"]["forceField"][subCls]["parameters"] = nonBondedList[subCls]["parameters"]
    simulation["topology"]["forceField"][subCls]["parameters"].update({"condition":"all"})
    simulation["topology"]["forceField"][subCls]["parameters"].update({"group":"pg{}".format(subCls)})
    if "labels" in nonBondedList[subCls].keys():
        simulation["topology"]["forceField"][subCls]["labels"] = nonBondedList[subCls]["labels"]
    if "data" in nonBondedList[subCls].keys():
        simulation["topology"]["forceField"][subCls]["data"] = nonBondedList[subCls]["data"]

    for i,inte in enumerate(integrators):
        simulation["simulationStep"][f"distances_{subCls}_{inte}"] = {}
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["type"] = ["ParticlesListMeasure", "DistancesMeasure"]
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["parameters"] = {}
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["parameters"]["startStep"] = i*nSteps
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["parameters"]["endStep"] = (i+1)*nSteps
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["parameters"]["intervalStep"] = nStepsOutput
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["parameters"]["outputFilePath"] = f"distances_{subCls}_{inte}.dat"
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["labels"] = ["id_i", "id_j"]
        simulation["simulationStep"][f"distances_{subCls}_{inte}"]["data"] = []

        for p in pairs:
            simulation["simulationStep"][f"distances_{subCls}_{inte}"]["data"].append(p)

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


