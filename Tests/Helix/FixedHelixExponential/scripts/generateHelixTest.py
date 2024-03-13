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

a        = parameters['a']
pitch    = parameters['pitch']
helicity = parameters['helicity']

sigma = parameters['sigma']

Kb = parameters['Kb']
Ka = parameters['Ka']
Kd = parameters['Kd']

E = parameters['E']

dt           = parameters['dt']
nSteps       = parameters['nSteps']
nStepsOutput = parameters['nStepsOutput']

angleRange = parameters.get('angleRange', 0.0)

print('Generating helix with parameters:')
print('N =', N)
print('a =', a)
print('pitch =', pitch)
print('helicity =', helicity)
print('sigma =', sigma)
print('Kb =', Kb)
print('Ka =', Ka)
print('Kd =', Kd)
print('E =', E)
print('angleRange =', angleRange)

R_H = computeHelixMatrix(a,pitch,helicity,sigma)
connectionNext,connectionPrevious = computeConnections(a,pitch,helicity,sigma)

print('R_H =\n', R_H)
print('connectionNext =', connectionNext)
print('connectionPrevious =', connectionPrevious)

pos,ori = generateHelix(N, a, pitch, helicity, sigma)
ori = helixAngularPerturbation(ori,angleRange)

#
simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "fixedHelixTest"

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
simulation["global"]["ensemble"]["data"]   = [[[1000.0,1000.0,1000.0], 1.0]]

simulation["integrator"] = {}

simulation["integrator"]["eulerMaruyamaRigid"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["type"] = ["Brownian", "EulerMaruyamaRigidBody"]
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["timeStep"]  = dt
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["viscosity"] = 1.0

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyamaRigid", nSteps],
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position","direction"]
simulation["state"]["data"] = []

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
simulation["topology"]["forceField"]["helix"]["type"]       = ["Bond2", "HelixExponential"]
simulation["topology"]["forceField"]["helix"]["parameters"] = {}

simulation["topology"]["forceField"]["helix"]["parameters"]["Kb"] = Kb
simulation["topology"]["forceField"]["helix"]["parameters"]["Ka"] = Ka
simulation["topology"]["forceField"]["helix"]["parameters"]["Kd"] = Kd

simulation["topology"]["forceField"]["helix"]["parameters"]["E"]  = E

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

simulation["topology"]["forceField"]["helix"]["parameters"]["e_x"] = e_x.tolist()
simulation["topology"]["forceField"]["helix"]["parameters"]["e_y"] = e_y.tolist()
simulation["topology"]["forceField"]["helix"]["parameters"]["e_z"] = e_z.tolist()

simulation["topology"]["forceField"]["helix"]["parameters"]["e_next"] = list(connectionNext)
simulation["topology"]["forceField"]["helix"]["parameters"]["e_prev"] = list(connectionPrevious)

simulation["topology"]["forceField"]["helix"]["labels"]     = ["id_i", "id_j"]
simulation["topology"]["forceField"]["helix"]["data"]       = []

for i in range(N-1):
    simulation["topology"]["forceField"]["helix"]["data"].append([i, i+1])

simulation["simulationStep"] = {}

simulation["simulationStep"]["pot"] = {}
simulation["simulationStep"]["pot"]["type"] = ["ParticlesListMeasure", "PotentialMeasure"]
simulation["simulationStep"]["pot"]["parameters"] = {}
simulation["simulationStep"]["pot"]["parameters"]["intervalStep"]   = nSteps
simulation["simulationStep"]["pot"]["parameters"]["outputFilePath"] = "pot.dat"
simulation["simulationStep"]["pot"]["labels"] = ["id"]
simulation["simulationStep"]["pot"]["data"]   = [[i] for i in range(N)]

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

simulation["simulationStep"]["output"] = {}
simulation["simulationStep"]["output"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["output"]["parameters"] = {}
simulation["simulationStep"]["output"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["output"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["output"]["parameters"]["outputFormat"]   = "sp"

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
