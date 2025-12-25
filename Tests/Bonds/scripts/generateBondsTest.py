import sys,os

import pyUAMMD

import numpy as np

import json
import jsbeautifier

from get_r0 import get_r0

#######################

COMPONENTS_PATH = os.path.join(os.getenv('UAMMD_PATH'),
                               'Components.json')

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

frictionConstant = param["frictionConstant"]
viscosity = param["viscosity"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

#Compute box size
L = 2.0*b*(N-1)
box = [L,L,L]

bondsList = param["bondsList"]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "bondsTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units",units]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

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

simulation["simulationStep"] = {}
simulation["simulationStep"]["groups"] = {}
simulation["simulationStep"]["groups"]["type"] = ["Groups", "GroupsList"]
simulation["simulationStep"]["groups"]["parameters"] = {}
simulation["simulationStep"]["groups"]["labels"] = ["name","type","selection"]
simulation["simulationStep"]["groups"]["data"] = []

offset = 0;
for simId,bnd in enumerate(components["Interactor"]["Bonds"]):

    cls,subCls,f = bnd

    if subCls not in bondsList:
        print("[WARNING] Bond type {} not found".format(subCls))
        continue

    simulation["simulationStep"]["groups"]["data"].append([f"pg{subCls}","BatchIds",[simId]])

    simulation["simulationStep"][f"write{subCls}"] = {}
    simulation["simulationStep"][f"write{subCls}"]["type"] = ["WriteStep", "WriteStep"]
    simulation["simulationStep"][f"write{subCls}"]["parameters"] = {"group" : f"pg{subCls}"}
    simulation["simulationStep"][f"write{subCls}"]["parameters"]["intervalStep"] = nStepsOutput
    simulation["simulationStep"][f"write{subCls}"]["parameters"]["outputFilePath"] = f"output{subCls}"
    simulation["simulationStep"][f"write{subCls}"]["parameters"]["outputFormat"] = "sp"

    r0 = get_r0(cls,subCls,bondsList)

    previous = np.asarray([r0,0,0])
    for i in range(N):
        simulation["state"]["data"].append([i+offset,list(previous)])
        #Generate random vector over a sphere of raidus r0
        theta = np.random.uniform(0,2*np.pi)
        phi   = np.random.uniform(0,np.pi)
        x = r0*np.sin(phi)*np.cos(theta)
        y = r0*np.sin(phi)*np.sin(theta)
        z = r0*np.cos(phi)

        previous = previous + np.asarray([x,y,z])


    for i in range(N):
        simulation["topology"]["structure"]["data"].append([i+offset, "A", simId])

    simulation["topology"]["forceField"][subCls] = {}
    simulation["topology"]["forceField"][subCls]["type"] = [cls, subCls]
    simulation["topology"]["forceField"][subCls]["parameters"] = bondsList[subCls]["commonParams"]

    bondParams   = bondsList[subCls]["bondParams"]

    needBounds   = bondsList[subCls]["needBounds"]

    simulation["topology"]["forceField"][subCls]["data"] = []
    if   cls == "Bond1":
        simulation["topology"]["forceField"][subCls]["labels"] = ["id_i"] + list(bondParams.keys())

        for i in range(N):
            d = [i+offset] + [bondParams[k] for k in bondParams.keys()]
            simulation["topology"]["forceField"][subCls]["data"].append(d)

    elif cls == "Bond2":
        simulation["topology"]["forceField"][subCls]["labels"] = ["id_i", "id_j"] + list(bondParams.keys())

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

        for i in range(N-1):
            d = [i+offset] + [i+offset+1]

            simulation["topology"]["forceField"][subCls]["data"].append(d + [bondParams[k] for k in bondParams.keys()])

            for i,inte in enumerate(integrators):
                simulation["simulationStep"][f"distances_{subCls}_{inte}"]["data"].append(d)

        if needBounds:

            boundParams = param["bond2bound"]

            tpy    = boundParams["type"]
            K      = boundParams["K"]
            factor = boundParams["factor"]

            simulation["topology"]["forceField"][subCls+"Bounds"] = {}
            simulation["topology"]["forceField"][subCls+"Bounds"]["type"] = ["Bond2", tpy]
            simulation["topology"]["forceField"][subCls+"Bounds"]["parameters"] = {}
            simulation["topology"]["forceField"][subCls+"Bounds"]["labels"] = ["id_i", "id_j","maxDistance","K"]
            simulation["topology"]["forceField"][subCls+"Bounds"]["data"] = []

            for i in range(N-1):
                d = [i+offset] + [i+offset+1]
                simulation["topology"]["forceField"][subCls+"Bounds"]["data"].append(d + [r0*factor,K])


    elif cls == "Bond3":
        simulation["topology"]["forceField"][subCls]["labels"] = ["id_i", "id_j", "id_k"] + list(bondParams.keys())

        for i,inte in enumerate(integrators):
            simulation["simulationStep"][f"angles_{subCls}_{inte}"] = {}
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["type"] = ["ParticlesListMeasure", "AnglesMeasure"]
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["parameters"] = {}
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["parameters"]["startStep"] = i*nSteps
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["parameters"]["endStep"] = (i+1)*nSteps
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["parameters"]["intervalStep"] = nStepsOutput
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["parameters"]["outputFilePath"] = f"angles_{subCls}_{inte}.dat"
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["labels"] = ["id_i", "id_j", "id_k"]
            simulation["simulationStep"][f"angles_{subCls}_{inte}"]["data"] = []

        for i in range(N-2):
            d = [i+offset] + [i+offset+1] + [i+offset+2]
            simulation["topology"]["forceField"][subCls]["data"].append(d + [bondParams[k] for k in bondParams.keys()])

            for i,inte in enumerate(integrators):
                simulation["simulationStep"][f"angles_{subCls}_{inte}"]["data"].append(d)

        simulation["topology"]["forceField"][subCls+"HarmonicBond2"] = {}
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["type"] = ["Bond2", "HarmonicCommon_K_r0"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["parameters"] = {"r0":bondsList["HarmonicCommon_K_r0"]["commonParams"]["r0"],
                                                                                      "K":bondsList["HarmonicCommon_K_r0"]["commonParams"]["K"]}
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["labels"] = ["id_i", "id_j"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["data"] = []

        for i in range(N-1):
            d = [i+offset] + [i+offset+1]
            simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["data"].append(d)

    elif cls == "Bond4":
        simulation["topology"]["forceField"][subCls]["labels"] = ["id_i", "id_j", "id_k", "id_l"] + list(bondParams.keys())

        for i,inte in enumerate(integrators):
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"] = {}
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["type"] = ["ParticlesListMeasure", "DihedralsMeasure"]
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["parameters"] = {}
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["parameters"]["startStep"] = i*nSteps
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["parameters"]["endStep"] = (i+1)*nSteps
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["parameters"]["intervalStep"] = nStepsOutput
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["parameters"]["outputFilePath"] = f"dihedrals_{subCls}_{inte}.dat"
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["labels"] = ["id_i", "id_j", "id_k", "id_l"]
            simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["data"] = []

        for i in range(N-3):
            d = [i+offset] + [i+offset+1] + [i+offset+2] + [i+offset+3]
            simulation["topology"]["forceField"][subCls]["data"].append(d + [bondParams[k] for k in bondParams.keys()])

            for i,inte in enumerate(integrators):
                simulation["simulationStep"][f"dihedrals_{subCls}_{inte}"]["data"].append(d)

        simulation["topology"]["forceField"][subCls+"HarmonicBond2"] = {}
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["type"] = ["Bond2", "HarmonicCommon_K_r0"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["parameters"] = {"r0":bondsList["HarmonicCommon_K_r0"]["commonParams"]["r0"],
                                                                                      "K":bondsList["HarmonicCommon_K_r0"]["commonParams"]["K"]}
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["labels"] = ["id_i", "id_j"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["data"] = []

        for i in range(N-1):
            d = [i+offset] + [i+offset+1]
            simulation["topology"]["forceField"][subCls+"HarmonicBond2"]["data"].append(d)

        simulation["topology"]["forceField"][subCls+"HarmonicBond3"] = {}
        simulation["topology"]["forceField"][subCls+"HarmonicBond3"]["type"] = ["Bond3", "HarmonicAngularCommon_K_theta0"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond3"]["parameters"] = {"theta0":bondsList["HarmonicAngularCommon_K_theta0"]["commonParams"]["theta0"],
                                                                                      "K":bondsList["HarmonicAngularCommon_K_theta0"]["commonParams"]["K"]}
        simulation["topology"]["forceField"][subCls+"HarmonicBond3"]["labels"] = ["id_i", "id_j", "id_k"]
        simulation["topology"]["forceField"][subCls+"HarmonicBond3"]["data"] = []

        for i in range(N-2):
            d = [i+offset] + [i+offset+1] + [i+offset+2]
            simulation["topology"]["forceField"][subCls+"HarmonicBond3"]["data"].append(d)
    else:
        print("Error: Unknown bond type: ", cls)
        sys.exit(1)

    offset += N

#Output

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
