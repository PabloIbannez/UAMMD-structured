import sys,os
import json

import numpy as np

import pyUAMMD

from scipy.spatial.transform import Rotation

with open('parameters.json', 'r') as f:
    parameters = json.load(f)

N = parameters['N']
c = parameters['c']

sigma = parameters['sigma']

E  = parameters['E']
K  = parameters['K']
rc = parameters['rc']

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']
nStepsPot    = parameters['nStepsPot']

avoidPBC     = parameters.get('avoidPBC', False)

connectionXup   = [ sigma/2.0,  0,  0]
connectionXdown = [-sigma/2.0,  0,  0]
connectionYup   = [ 0,  sigma/2.0,  0]
connectionYdown = [ 0, -sigma/2.0,  0]
connectionZup   = [ 0,  0,  sigma/2.0]
connectionZdown = [ 0,  0, -sigma/2.0]

print('Generating cosine with parameters:')
print('N =', N)
print('c =', c)
print('sigma =', sigma)
print('E =', E)
print('K =', K)
print('rc =', rc)

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
simulation["system"]["info"]["parameters"]["name"] = "dynamicCosineTest"

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

for i in range(N):
    # Generate random position within the box
    x = np.random.uniform(-L/2,L/2)
    y = np.random.uniform(-L/2,L/2)
    z = np.random.uniform(-L/2,L/2)
    p = [x,y,z]

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

simulation["topology"]["forceField"]["cosine"] = {}
simulation["topology"]["forceField"]["cosine"]["type"]       = ["PatchyParticles", "DynamicallyBondedPatchyParticles"]

simulation["topology"]["forceField"]["cosine"]["patchesState"] = {}
simulation["topology"]["forceField"]["cosine"]["patchesState"]["labels"] = ["id", "position"]
simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"]   = []

for i in range(N):
    index=i*6
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index  ), list(connectionXup)])
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index+1), list(connectionXdown)])
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index+2), list(connectionYup)])
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index+3), list(connectionYdown)])
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index+4), list(connectionZup)])
    simulation["topology"]["forceField"]["cosine"]["patchesState"]["data"].append([int(index+5), list(connectionZdown)])

simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]={}
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["fundamental"] = {}
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["fundamental"]["type"]       = ["Fundamental","DynamicallyBondedPatchyParticles"]
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["fundamental"]["parameters"] = {"energyThreshold":0.0}

simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["types"]  = {}
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["types"]["type"]   = ["Types","Basic"]
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["topology"]["forceField"]["cosine"]["patchesGlobal"]["types"]["data"]   = [["P",0.0,sigma/10.0,0.0]]

simulation["topology"]["forceField"]["cosine"]["patchesTopology"]={}

simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]={}
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["labels"] = ["id", "type", "parentId"]
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"]   = []

for i in range(N):
    index=i*6
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index  ,"P",i])
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index+1,"P",i])
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index+2,"P",i])
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index+3,"P",i])
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index+4,"P",i])
    simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["structure"]["data"].append([index+5,"P",i])

simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"] = {}

#Verlet list
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["verletList"]={}
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["verletList"]["type"]       =  ["VerletConditionalListSet", "all"]
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["verletList"]["parameters"] =  {"cutOffVerletFactor":3.0}

simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]={}
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]["type"]       =  ["NonBondedPatches", "DistanceSwitchCosine"]
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]["parameters"] =  {"condition":"all"}
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]["labels"]     =  ["name_i", "name_j", "E", "K", "rc"]
simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]["data"]       =  []

d = ["P", "P", E, K, rc]

simulation["topology"]["forceField"]["cosine"]["patchesTopology"]["forceField"]["cosine"]["data"].append(d)

#######################################################

simulation["simulationStep"] = {}

simulation["simulationStep"]["pot"] = {}
simulation["simulationStep"]["pot"]["type"] = ["ParticlesListMeasure", "PotentialMeasure"]
simulation["simulationStep"]["pot"]["parameters"] = {}
simulation["simulationStep"]["pot"]["parameters"]["intervalStep"]   = nStepsPot
simulation["simulationStep"]["pot"]["parameters"]["outputFilePath"] = "pot.dat"
simulation["simulationStep"]["pot"]["labels"] = ["id"]
simulation["simulationStep"]["pot"]["data"]   = [[i] for i in range(N)]

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
