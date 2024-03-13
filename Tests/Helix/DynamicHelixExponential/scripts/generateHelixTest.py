import sys,os
import json

import numpy as np

import pyUAMMD

# We want to import functions defined at ../scripts/helixGeneration.py
sys.path.append(os.path.join(sys.path[0],'..','..','scripts'))
from helixGeneration import computeHelixMatrix,computeConnections,generateHelix,helixAngularPerturbation

from scipy.spatial.transform import Rotation

with open('parameters.json', 'r') as f:
    parameters = json.load(f)

N = parameters['N']
c = parameters['c']

a        = parameters['a']
pitch    = parameters['pitch']
helicity = parameters['helicity']

sigma = parameters['sigma']

Kb = parameters['Kb']
Ka = parameters['Ka']
Kd = parameters['Kd']

E = parameters['E']

rc = parameters['rc']

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']
nStepsPot    = parameters['nStepsPot']

initHelix  = parameters.get('initHelix', False)
angleRange = parameters.get('angleRange', 0.0)

avoidPBC   = parameters.get('avoidPBC', False)

if not initHelix:
    if "angleRange" in parameters:
        print("WARNING: angleRange is ignored when initHelix is False")

print('Generating helix with parameters:')
print('N =', N)
print('c =', c)
print('a =', a)
print('pitch =', pitch)
print('helicity =', helicity)
print('sigma =', sigma)
print('Kb =', Kb)
print('Ka =', Ka)
print('Kd =', Kd)
print('E =', E)
print('angleRange =', angleRange)

# c = N/L^3
L = (N/c)**(1.0/3.0)
box = [L,L,L]

if avoidPBC:
    box[0] = 2.0*L
    box[1] = 2.0*L
    box[2] = 2.0*L

R_H = computeHelixMatrix(a,pitch,helicity,sigma)
connectionNext,connectionPrevious = computeConnections(a,pitch,helicity,sigma)

pos,ori = generateHelix(N, a, pitch, helicity, sigma)
ori = helixAngularPerturbation(ori,angleRange)

print('R_H =\n', R_H)
print('connectionNext =', connectionNext)
print('connectionPrevious =', connectionPrevious)

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "dynamicHelixTest"

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

if not initHelix:
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
else:
    for i in range(N):

        p = list(pos[i])

        o = ori[i]
        q = Rotation.from_matrix(o).as_quat() # scalar last
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

simulation["topology"]["forceField"]["helix"] = {}
simulation["topology"]["forceField"]["helix"]["type"]       = ["PatchyParticles", "DynamicallyBondedPatchyParticles"]

simulation["topology"]["forceField"]["helix"]["patchesState"] = {}
simulation["topology"]["forceField"]["helix"]["patchesState"]["labels"] = ["id", "position"]
simulation["topology"]["forceField"]["helix"]["patchesState"]["data"]   = []

for i in range(N):
    index=i*2
    simulation["topology"]["forceField"]["helix"]["patchesState"]["data"].append([int(index  ), list(connectionPrevious)])
    simulation["topology"]["forceField"]["helix"]["patchesState"]["data"].append([int(index+1), list(connectionNext)])

simulation["topology"]["forceField"]["helix"]["patchesGlobal"]={}
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["fundamental"] = {}
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["fundamental"]["type"]       = ["Fundamental","DynamicallyBondedPatchyParticles"]
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["fundamental"]["parameters"] = {"energyThreshold":0.0}

simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["types"]  = {}
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["types"]["type"]   = ["Types","Basic"]
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["topology"]["forceField"]["helix"]["patchesGlobal"]["types"]["data"]   = [["S",0.0,sigma/10.0,0.0],
                                                                                     ["E",0.0,sigma/10.0,0.0]]

simulation["topology"]["forceField"]["helix"]["patchesTopology"]={}

simulation["topology"]["forceField"]["helix"]["patchesTopology"]["structure"]={}
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["structure"]["labels"] = ["id", "type", "parentId"]
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["structure"]["data"]   = []

for i in range(N):
    index=i*2
    simulation["topology"]["forceField"]["helix"]["patchesTopology"]["structure"]["data"].append([index  ,"E",i])
    simulation["topology"]["forceField"]["helix"]["patchesTopology"]["structure"]["data"].append([index+1,"S",i])

simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"] = {}

#Verlet list
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["verletList"]={}
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["verletList"]["type"]       =  ["VerletConditionalListSet", "interDifferentType"]
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["verletList"]["parameters"] =  {"cutOffVerletFactor":3.0}

simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]={}
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["type"]       =  ["NonBondedPatches", "HelixExponential"]
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["parameters"] =  {"condition":"inter"}
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["labels"]     =  ["name_i", "name_j", "Eb", "Kb", "Ka", "Kd", "rc", "e_x", "e_y", "e_z"]
simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["data"]       =  []

d = ["E", "S", E, Kb, Ka, Kd, rc]

e_x = R_H[:,0]
e_y = R_H[:,1]
e_z = R_H[:,2]

e_x_norm = np.linalg.norm(e_x)
e_y_norm = np.linalg.norm(e_y)
e_z_norm = np.linalg.norm(e_z)

if np.abs(e_x_norm - 1.0) > 1e-6:
    print('e_x is not normalized!')
    print('e_x_norm =', e_x_norm)
    exit()
if np.abs(e_y_norm - 1.0) > 1e-6:
    print('e_y is not normalized!')
    print('e_y_norm =', e_y_norm)
    exit()
if np.abs(e_z_norm - 1.0) > 1e-6:
    print('e_z is not normalized!')
    print('e_z_norm =', e_z_norm)
    exit()

e_xy = np.dot(e_x,e_y)
e_xz = np.dot(e_x,e_z)
e_yz = np.dot(e_y,e_z)

if abs(e_xy) > 1e-6:
    print('e_x and e_y are not orthogonal!')
    print('e_xy =', e_xy)
    sys.exit()

if abs(e_xz) > 1e-6:
    print('e_x and e_z are not orthogonal!')
    print('e_xz =', e_xz)
    sys.exit()

if abs(e_yz) > 1e-6:
    print('e_y and e_z are not orthogonal!')
    print('e_yz =', e_yz)
    sys.exit()

d.append(e_x.tolist())
d.append(e_y.tolist())
d.append(e_z.tolist())

simulation["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["data"].append(d)

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
