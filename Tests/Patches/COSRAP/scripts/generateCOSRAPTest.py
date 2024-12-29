import sys,os
import json

import numpy as np
import quaternion as qt

import pyUAMMD

from scipy.spatial.transform import Rotation

with open('parameters.json', 'r') as f:
    parameters = json.load(f)

N = parameters['N']
c = parameters['c']

sigma = parameters['sigma']

E    = parameters['E']
Kswt = parameters['Kswt']
Krap = parameters['Krap']
rc   = parameters['rc']

Rq    = parameters['R']
Rqinv = qt.quaternion(*Rq).conjugate()
Rqinv = [Rqinv.w, Rqinv.x, Rqinv.y, Rqinv.z]

## Convert Rq to rotation matrix
#R    = qt.as_rotation_matrix(qt.quaternion(*Rq))
#Rinv = qt.as_rotation_matrix(qt.quaternion(*Rqinv))
#
## Check
#I = R@Rinv
#print('I =', I)

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']
nStepsPot    = parameters['nStepsPot']

avoidPBC     = parameters.get('avoidPBC', False)

connectionS = [ sigma/2.0,  0,  0]
connectionE = [-sigma/2.0,  0,  0]

print('Generating COSRAP with parameters:')
print('N =', N)
print('c =', c)
print('sigma =', sigma)
print('E =', E)
print('Kswt =', Kswt)
print('Krap =', Krap)
print('Rq =', Rq)
print('Rqinv =', Rqinv)
print('rc =', rc)

# c = N/L^3
L = (N/c)**(1.0/3.0)
box = [L,L,L]

if avoidPBC:
    box[0] = L+5.0*sigma
    box[1] = L+5.0*sigma
    box[2] = L+5.0*sigma

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "COSRAPTest"

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
simulation["integrator"]["schedule"]["type"]   = ["Schedule", "Integrator"]
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

simulation["topology"]["forceField"]["cosrap"] = {}
simulation["topology"]["forceField"]["cosrap"]["type"]       = ["PatchyParticles", "DynamicallyBondedPatchyParticles"]

simulation["topology"]["forceField"]["cosrap"]["patchesState"] = {}
simulation["topology"]["forceField"]["cosrap"]["patchesState"]["labels"] = ["id", "position"]
simulation["topology"]["forceField"]["cosrap"]["patchesState"]["data"]   = []

for i in range(N):
    index=i*2
    simulation["topology"]["forceField"]["cosrap"]["patchesState"]["data"].append([int(index  ), list(connectionS)])
    simulation["topology"]["forceField"]["cosrap"]["patchesState"]["data"].append([int(index+1), list(connectionE)])

simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]={}
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["fundamental"] = {}
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["fundamental"]["type"]       = ["Fundamental","DynamicallyBondedPatchyParticles"]
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["fundamental"]["parameters"] = {"energyThreshold":0.0}

simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["types"]  = {}
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["types"]["type"]   = ["Types","Basic"]
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["topology"]["forceField"]["cosrap"]["patchesGlobal"]["types"]["data"]   = [["S",0.0,sigma/10.0,0.0],
                                                                                      ["E",0.0,sigma/10.0,0.0]]

simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]={}

simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["structure"]={}
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["structure"]["labels"] = ["id", "type", "parentId"]
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["structure"]["data"]   = []

for i in range(N):
    index=i*2
    simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["structure"]["data"].append([index  ,"S",i])
    simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["structure"]["data"].append([index+1,"E",i])

simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"] = {}

#Verlet list
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["verletList"]={}
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["verletList"]["type"]       =  ["VerletConditionalListSet", "interDifferentType"]
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["verletList"]["parameters"] =  {"cutOffVerletFactor":1.1}

simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]={}
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["type"]       =  ["NonBondedPatches", "COSRAP"]
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["parameters"] =  {"condition":"inter"}
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["labels"]     =  ["name_i", "name_j", "E", "Kswt", "Krap", "rc", "R"]
simulation["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["data"]       =  [["S","E",E,Kswt,Krap,rc,Rq],
                                                                                                            ["E","S",E,Kswt,Krap,rc,Rqinv]]

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
