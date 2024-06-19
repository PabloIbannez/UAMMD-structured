import sys,os
import json

from tqdm import tqdm

import numpy as np

import pyUAMMD

from scipy.spatial.transform import Rotation

with open('parameters.json', 'r') as f:
    parameters = json.load(f)

N = parameters['N']
c = parameters['c']

sigma = parameters['sigma']

E  = parameters['E']
rc = parameters['rc']
Kb = parameters['Kb']

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']

avoidPBC     = parameters.get('avoidPBC', False)

connectionXup   = [ sigma/2.0,  0,  0]
connectionXdown = [-sigma/2.0,  0,  0]
connectionYup   = [ 0,  sigma/2.0,  0]
connectionYdown = [ 0, -sigma/2.0,  0]
connectionZup   = [ 0,  0,  sigma/2.0]
connectionZdown = [ 0,  0, -sigma/2.0]

print('Generating expo with parameters:')
print('N =', N)
print('c =', c)
print('sigma =', sigma)
print('E =', E)
print('rc =', rc)
print('Kb =', Kb)

# c = N/L^3
L = (N/c)**(1.0/3.0)
box = [L,L,L]

if avoidPBC:
    box[0] = 2.0*L
    box[1] = 2.0*L
    box[2] = 2.0*L

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "dynamicExponentialTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","None"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", 1.0, sigma/2, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, 1.0]]

simulation["integrator"] = {}

if dt == 0.0:
    simulation["integrator"]["integrator"] = {}
    simulation["integrator"]["integrator"]["type"] = ["None", "None"]
else:
    simulation["integrator"]["integrator"] = {}
    simulation["integrator"]["integrator"]["type"] = ["Brownian", "EulerMaruyamaRigidBody"]
    simulation["integrator"]["integrator"]["parameters"] = {}
    simulation["integrator"]["integrator"]["parameters"]["timeStep"]  = dt
    simulation["integrator"]["integrator"]["parameters"]["viscosity"] = 1.0


simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "integrator", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position","direction"]
simulation["state"]["data"] = []

for i in tqdm(range(N)):
    # Generate random position within the box
    clash = True
    while clash:
        x = np.random.uniform(-L/2,L/2)
        y = np.random.uniform(-L/2,L/2)
        z = np.random.uniform(-L/2,L/2)
        p = [x,y,z]

        if i == 0:
            clash = False

        for j in range(i):
            dx = p[0] - simulation["state"]["data"][j][1][0]
            dy = p[1] - simulation["state"]["data"][j][1][1]
            dz = p[2] - simulation["state"]["data"][j][1][2]

            #Consider PBC
            dx = dx - box[0]*np.round(dx/box[0])
            dy = dy - box[1]*np.round(dy/box[1])
            dz = dz - box[2]*np.round(dz/box[2])

            r2 = dx*dx + dy*dy + dz*dz

            if r2 < 0.8*sigma*sigma:
                clash = True
                break
            else:
                clash = False

    q = Rotation.random().as_quat() # scalar last
    q = [q[3], q[0], q[1], q[2]] # scalar first
    q = list(q)

    simulation["state"]["data"].append([i, p, q])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"]   = []

for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

simulation["topology"]["forceField"] = {}

#Verlet list
simulation["topology"]["forceField"]["nl"]={}
simulation["topology"]["forceField"]["nl"]["type"]       =  ["VerletConditionalListSet", "all"]
simulation["topology"]["forceField"]["nl"]["parameters"] =  {"cutOffVerletFactor":1.2}

simulation["topology"]["forceField"]["wca"] = {}
simulation["topology"]["forceField"]["wca"]["type"]       = ["NonBonded", "WCAType2"]
simulation["topology"]["forceField"]["wca"]["parameters"] = {"cutOffFactor":2.5,
                                                             "condition":"all"}
simulation["topology"]["forceField"]["wca"]["labels"]     = ["name_i", "name_j", "epsilon", "sigma"]
simulation["topology"]["forceField"]["wca"]["data"]       = [["A", "A", 1.0, sigma]]

simulation["topology"]["forceField"]["exponential"] = {}
simulation["topology"]["forceField"]["exponential"]["type"]       = ["PatchyParticles", "DynamicallyBondedPatchyParticles"]

simulation["topology"]["forceField"]["exponential"]["patchesState"] = {}
simulation["topology"]["forceField"]["exponential"]["patchesState"]["labels"] = ["id", "position"]
simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"]   = []

for i in range(N):
    index=i*6
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index  ), list(connectionXup)])
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index+1), list(connectionXdown)])
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index+2), list(connectionYup)])
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index+3), list(connectionYdown)])
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index+4), list(connectionZup)])
    simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([int(index+5), list(connectionZdown)])

simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]={}
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["fundamental"] = {}
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["fundamental"]["type"]       = ["Fundamental","DynamicallyBondedPatchyParticles"]
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["fundamental"]["parameters"] = {"energyThreshold":0.0}

simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["types"]  = {}
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["types"]["type"]   = ["Types","Basic"]
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["topology"]["forceField"]["exponential"]["patchesGlobal"]["types"]["data"]   = [["P",0.0,sigma/10.0,0.0]]

simulation["topology"]["forceField"]["exponential"]["patchesTopology"]={}

simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]={}
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["labels"] = ["id", "type", "parentId"]
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"]   = []

for i in range(N):
    index=i*6
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index  ,"P",i])
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index+1,"P",i])
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index+2,"P",i])
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index+3,"P",i])
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index+4,"P",i])
    simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([index+5,"P",i])

simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"] = {}

#Verlet list
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["verletList"]={}
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["verletList"]["type"]       =  ["VerletConditionalListSet", "all"]
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["verletList"]["parameters"] =  {"cutOffVerletFactor":1.2}

simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]={}
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]["type"]       =  ["NonBondedPatches", "DistanceSwitchExponential"]
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]["parameters"] =  {"condition":"all"}
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]["labels"]     =  ["name_i", "name_j", "E", "K", "rc"]
simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]["data"]       =  []

d = ["P", "P", E, Kb, rc]

simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["forceField"]["exponential"]["data"].append(d)

#######################################################

simulation["simulationStep"] = {}

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

simulation["simulationStep"]["output"] = {}
simulation["simulationStep"]["output"]["type"] = ["WriteStep", "WritePatchyParticlesStep"]
simulation["simulationStep"]["output"]["parameters"] = {}
simulation["simulationStep"]["output"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"]["output"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["output"]["parameters"]["outputFormat"]   = "sp"

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
