import os
import json

import numpy as np

import icosphere

from pyUAMMD import simulation

with open('parameters.json') as f:
    parameters = json.load(f)

partRadius    = parameters['particleRadius']
resolution    = parameters['resolution']
R             = parameters['radius']

Kb            = parameters['Kb']
Kf            = parameters['Kf']

L             = parameters['boxSize']
temperature   = parameters['temperature']

timeStep      = parameters['timeStep']
friction      = parameters['frictionConstant']

lambd         = parameters['lambda']
n             = parameters['n']
alpha         = parameters['alpha']

epsilon       = parameters['epsilon']
sigma         = parameters['sigma']

nSteps        = 1
nStepsOutput  = 1
nStepsMeasure = 1

vertices, faces = icosphere.icosphere(nu = resolution)

ids = list(range(vertices.shape[0]))

for i in range(len(vertices)):
    vertices[i] = vertices[i]*R

dst = np.linalg.norm(vertices[1] - vertices[0])

#For each vertex sum small random displacement
for i in range(len(vertices)):
    vertices[i] = vertices[i] + np.random.rand(3)*partRadius*1.5

fixed = []
for v in vertices:
    vf = v+np.random.rand(3)*1.05
    fixed.append([vf[0], vf[1], vf[2]])

bonds = set()
for face in faces:
    bonds.add(tuple(sorted([face[0], face[1]])))
    bonds.add(tuple(sorted([face[1], face[2]])))
    bonds.add(tuple(sorted([face[2], face[0]])))

sim = simulation(debug = True)

sim["system"] = {}
sim["system"]["info"] = {}
sim["system"]["info"]["type"] = ["Simulation","Information"]
sim["system"]["info"]["parameters"] = {}
sim["system"]["info"]["parameters"]["name"] = "tiTest"

sim["global"] = {}

sim["global"]["units"] = {}
sim["global"]["units"]["type"] = ["Units","None"]

sim["global"]["fundamental"] = {}
sim["global"]["fundamental"]["type"] = ["Fundamental","Time"]
sim["global"]["fundamental"]["parameters"] = {}
sim["global"]["fundamental"]["parameters"]["currentStep"] = 0
sim["global"]["fundamental"]["parameters"]["simulationTime"] = 0.0

sim["global"]["types"] = {}
sim["global"]["types"]["type"]   = ["Types","Basic"]
sim["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
sim["global"]["types"]["data"]   = [["A", 1.0, partRadius, 0.0]]

sim["global"]["ensemble"] = {}
sim["global"]["ensemble"]["type"]   = ["Ensemble","NVTlambda"]
sim["global"]["ensemble"]["labels"] = ["box", "temperature", "lambda"]
sim["global"]["ensemble"]["data"]   = [[[L, L, L], temperature,lambd]]

sim["integrator"] = {}

sim["integrator"]["bbk"] = {}
sim["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
sim["integrator"]["bbk"]["parameters"] = {}
sim["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
sim["integrator"]["bbk"]["parameters"]["frictionConstant"] = friction

sim["integrator"]["schedule"] = {}
sim["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
sim["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
sim["integrator"]["schedule"]["data"] = [
    [1, "bbk", nSteps]
]

sim["state"] = {}
sim["state"]["labels"] = ["id", "position"]
sim["state"]["data"] = []
for i in range(vertices.shape[0]):
    pos = vertices[i]
    sim["state"]["data"].append([i, list(pos)])

sim["topology"] = {}
sim["topology"]["structure"] = {}
sim["topology"]["structure"]["labels"] = ["id", "type","batchId"]
sim["topology"]["structure"]["data"] = []
for i in range(vertices.shape[0]):
    sim["topology"]["structure"]["data"].append([i, "A", 0])

sim["topology"]["forceField"] = {}

sim["topology"]["forceField"]["fixed"] = {}
sim["topology"]["forceField"]["fixed"]["type"]       = ["Bond1","LambdaFixedHarmonic"]
sim["topology"]["forceField"]["fixed"]["labels"]     = ["id_i","r0","K","position"]
sim["topology"]["forceField"]["fixed"]["parameters"] = {"n":n}
sim["topology"]["forceField"]["fixed"]["data"]       = []

for i,fi in enumerate(fixed):
    sim["topology"]["forceField"]["fixed"]["data"].append([i,[0.0,0.0,0.0], [Kf,Kf,Kf], [fi[0], fi[1], fi[2]]])

sim["topology"]["forceField"]["bonds"] = {}
sim["topology"]["forceField"]["bonds"]["type"]       = ["Bond2","LambdaHarmonic"]
sim["topology"]["forceField"]["bonds"]["labels"]     = ["id_i","id_j","r0","K"]
sim["topology"]["forceField"]["bonds"]["parameters"] = {"n":n}
sim["topology"]["forceField"]["bonds"]["data"]       = []

for bond in bonds:
    p1  = vertices[bond[0]]
    p2  = vertices[bond[1]]
    sim["topology"]["forceField"]["bonds"]["data"].append([int(bond[0]), int(bond[1]), dst, Kb])

sim["topology"]["forceField"]["ljbonds"] = {}
sim["topology"]["forceField"]["ljbonds"]["type"]       = ["Bond2","LennardJonesSoftCoreType2"]
sim["topology"]["forceField"]["ljbonds"]["labels"]     = ["id_i","id_j","epsilon","sigma"]
sim["topology"]["forceField"]["ljbonds"]["parameters"] = {"alpha":alpha, "n":n}
sim["topology"]["forceField"]["ljbonds"]["data"]       = []

for bond in bonds:
    p1  = vertices[bond[0]]
    p2  = vertices[bond[1]]
    sim["topology"]["forceField"]["ljbonds"]["data"].append([int(bond[0]), int(bond[1]), epsilon, dst])

sim["topology"]["forceField"]["nl"] = {}
sim["topology"]["forceField"]["nl"]["type"] = ["VerletConditionalListSet","all"]
sim["topology"]["forceField"]["nl"]["parameters"] = {"cutOffVerletFactor": 1.2}

sim["topology"]["forceField"]["lj"] = {}
sim["topology"]["forceField"]["lj"]["type"]       = ["NonBonded","LennardJonesSoftCoreType2"]
sim["topology"]["forceField"]["lj"]["parameters"] = {"cutOffFactor":2.5,"condition":"all","alpha":alpha, "n":n}
sim["topology"]["forceField"]["lj"]["labels"]     = ["name_i","name_j","epsilon","sigma"]
sim["topology"]["forceField"]["lj"]["data"]       = [["A","A",epsilon,sigma]]

#Output

sim["simulationStep"] = {}

sim["simulationStep"]["info"] = {}
sim["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
sim["simulationStep"]["info"]["parameters"] = {}
sim["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

sim["simulationStep"]["write"] = {}
sim["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
sim["simulationStep"]["write"]["parameters"] = {}
sim["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
sim["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
sim["simulationStep"]["write"]["parameters"]["outputFormat"] = "sp"

sim["simulationStep"]["pot"] = {}
sim["simulationStep"]["pot"]["type"] = ["ParticlesListMeasure","PotentialMeasure"]
sim["simulationStep"]["pot"]["parameters"] = {}
sim["simulationStep"]["pot"]["parameters"]["intervalStep"] = nStepsOutput
sim["simulationStep"]["pot"]["parameters"]["outputFilePath"] = "pot.dat"
sim["simulationStep"]["pot"]["labels"] = ["id"]
sim["simulationStep"]["pot"]["data"] = [[i] for i in ids]

sim["simulationStep"]["ti"] = {}
sim["simulationStep"]["ti"]["type"] = ["ThermodynamicMeasure","ThermodynamicIntegration"]
sim["simulationStep"]["ti"]["parameters"] = {}
sim["simulationStep"]["ti"]["parameters"]["intervalStep"] = nStepsOutput
sim["simulationStep"]["ti"]["parameters"]["stepLambda"]   = 1
sim["simulationStep"]["ti"]["parameters"]["initLambda"]   = lambd
sim["simulationStep"]["ti"]["parameters"]["finalLambda"]  = 0.0
sim["simulationStep"]["ti"]["parameters"]["lambdaIntervalLength"] = 10
sim["simulationStep"]["ti"]["parameters"]["outputFilePath"] = "ti.dat"

#Create results folder if it does not exist
if not os.path.exists("results"):
    os.makedirs("results")

sim.write("results/simulation.json")



