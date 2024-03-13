import sys,os
import json

import numpy as np

import pyUAMMD

from scipy.spatial.transform import Rotation

with open('parameters.json', 'r') as f:
    parameters = json.load(f)

N = parameters['N']
c = parameters['c']

E = parameters['E']
D = parameters['D']

linkerPos = parameters['linkerPos']

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']
nStepsPot    = parameters['nStepsPot']

print('Generating linker particles')
print('N =', N)
print('c =', c)
print('E =', E)
print('D =', D)
print('linkerPos =', linkerPos)

# c = N/L^3
L = (N/c)**(1.0/3.0)
box = [L,L,L]

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "linkerTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","None"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", 1.0, 0.5, 0.0]]

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
    x = np.random.uniform(-0.5,0.5)
    y = np.random.uniform(-0.5,0.5)
    z = np.random.uniform(-0.5,0.5)
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

simulation["topology"]["forceField"]["linker"] = {}
simulation["topology"]["forceField"]["linker"]["type"]       = ["PatchyParticles", "PatchyParticles"]

simulation["topology"]["forceField"]["linker"]["patchesState"] = {}
simulation["topology"]["forceField"]["linker"]["patchesState"]["labels"] = ["id", "position"]
simulation["topology"]["forceField"]["linker"]["patchesState"]["data"]   = []

for i in range(N):
    simulation["topology"]["forceField"]["linker"]["patchesState"]["data"].append([i,linkerPos])

simulation["topology"]["forceField"]["linker"]["patchesGlobal"]={}
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["fundamental"] = {}
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["fundamental"]["type"]       = ["Fundamental","DynamicallyBondedPatchyParticles"]
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["fundamental"]["parameters"] = {"energyThreshold":0.0}

simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["types"]  = {}
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["types"]["type"]   = ["Types","Basic"]
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["topology"]["forceField"]["linker"]["patchesGlobal"]["types"]["data"]   = [["L",0.0,1.0/10.0,0.0]]

simulation["topology"]["forceField"]["linker"]["patchesTopology"]={}

simulation["topology"]["forceField"]["linker"]["patchesTopology"]["structure"]={}
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["structure"]["labels"] = ["id", "type", "parentId"]
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["structure"]["data"]   = []

for i in range(N):
    simulation["topology"]["forceField"]["linker"]["patchesTopology"]["structure"]["data"].append([i,"L",i])

simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"] = {}

simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]={}
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["type"]       =  ["SurfacePatches", "Linker"]
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["parameters"] =  {}
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["labels"]     =  ["name", "epsilon", "sigma"]
simulation["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["data"]       =  [["L", E, D]]

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
