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
Kd            = parameters['Kd']
Kf            = parameters['Kf']
Kc            = parameters['Kc']

L             = parameters['boxSize']
temperature   = parameters['temperature']

timeStep      = parameters['timeStep']
friction      = parameters['frictionConstant']

nSteps        = parameters['nSteps']
nStepsOutput  = parameters['nStepsOutput']
nStepsMeasure = parameters['nStepsMeasure']

separation = 2*R*1.5

vertices1, faces1 = icosphere.icosphere(nu = resolution)
vertices2, faces2 = icosphere.icosphere(nu = resolution)
vertices3, faces3 = icosphere.icosphere(nu = resolution)

idOffset2 = vertices1.shape[0]
idOffset3 = vertices1.shape[0] + vertices2.shape[0]

ids1 = list(range(vertices1.shape[0]))
ids2 = list(range(idOffset2, idOffset2 + vertices2.shape[0]))
ids3 = list(range(idOffset3, idOffset3 + vertices3.shape[0]))

for i in range(len(vertices1)):
    vertices1[i] = vertices1[i]*R + np.array([0, 0, L/2])

for i in range(len(vertices2)):
    vertices2[i] = vertices2[i]*R

for i in range(len(vertices3)):
    vertices3[i] = vertices3[i]*R + np.array([separation, 0, 0])

bonds1 = set()
for face in faces1:
    bonds1.add(tuple(sorted([face[0], face[1]])))
    bonds1.add(tuple(sorted([face[1], face[2]])))
    bonds1.add(tuple(sorted([face[2], face[0]])))

bonds2 = set()
for face in faces2:
    bonds2.add(tuple(sorted([face[0] + idOffset2, face[1] + idOffset2])))
    bonds2.add(tuple(sorted([face[1] + idOffset2, face[2] + idOffset2])))
    bonds2.add(tuple(sorted([face[2] + idOffset2, face[0] + idOffset2])))

bonds3 = set()
for face in faces3:
    bonds3.add(tuple(sorted([face[0] + idOffset3, face[1] + idOffset3])))
    bonds3.add(tuple(sorted([face[1] + idOffset3, face[2] + idOffset3])))
    bonds3.add(tuple(sorted([face[2] + idOffset3, face[0] + idOffset3])))

#Dihedral angles
facesPairs1 = set()
for i in range(len(faces1)):
    for j in range(i+1, len(faces1)):
        if len(set(faces1[i]).intersection(set(faces1[j]))) == 2:
            facesPairs1.add(tuple(sorted([i, j])))

facesPairs2 = set()
for i in range(len(faces2)):
    for j in range(i+1, len(faces2)):
        if len(set(faces2[i]).intersection(set(faces2[j]))) == 2:
            facesPairs2.add(tuple(sorted([i, j])))

facesPairs3 = set()
for i in range(len(faces3)):
    for j in range(i+1, len(faces3)):
        if len(set(faces3[i]).intersection(set(faces3[j]))) == 2:
            facesPairs3.add(tuple(sorted([i, j])))

dihedrals1 = []
for pair in facesPairs1:
    face1 = faces1[pair[0]]
    face2 = faces1[pair[1]]

    common = set(face1).intersection(set(face2))

    if len(common) == 2:
        face1Id = list(set(face1).difference(common))[0]
        face2Id = list(set(face2).difference(common))[0]

        i = int(face1Id)
        j = int(common.pop())
        k = int(common.pop())
        l = int(face2Id)

        dihedrals1.append([i, j, k, l])

dihedrals2 = []
for pair in facesPairs2:
    face1 = faces2[pair[0]]
    face2 = faces2[pair[1]]

    common = set(face1).intersection(set(face2))

    if len(common) == 2:
        face1Id = list(set(face1).difference(common))[0]
        face2Id = list(set(face2).difference(common))[0]

        i = int(face1Id) + idOffset2
        j = int(common.pop()) + idOffset2
        k = int(common.pop()) + idOffset2
        l = int(face2Id) + idOffset2

        dihedrals2.append([i, j, k, l])

dihedrals3 = []
for pair in facesPairs3:
    face1 = faces3[pair[0]]
    face2 = faces3[pair[1]]

    common = set(face1).intersection(set(face2))

    if len(common) == 2:
        face1Id = list(set(face1).difference(common))[0]
        face2Id = list(set(face2).difference(common))[0]

        i = int(face1Id) + idOffset3
        j = int(common.pop()) + idOffset3
        k = int(common.pop()) + idOffset3
        l = int(face2Id) + idOffset3

        dihedrals3.append([i, j, k, l])

sim = simulation(debug = True)

sim["system"] = {}
sim["system"]["info"] = {}
sim["system"]["info"]["type"] = ["Simulation","Information"]
sim["system"]["info"]["parameters"] = {}
sim["system"]["info"]["parameters"]["name"] = "setsTest"

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
sim["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
sim["global"]["ensemble"]["labels"] = ["box", "temperature"]
sim["global"]["ensemble"]["data"]   = [[[L, L, L], temperature]]

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
for i in range(vertices1.shape[0]):
    pos = vertices1[i]
    sim["state"]["data"].append([i, list(pos)])
for i in range(vertices2.shape[0]):
    pos = vertices2[i]
    sim["state"]["data"].append([i + idOffset2, list(pos)])
for i in range(vertices3.shape[0]):
    pos = vertices3[i]
    sim["state"]["data"].append([i + idOffset3, list(pos)])

sim["topology"] = {}
sim["topology"]["structure"] = {}
sim["topology"]["structure"]["labels"] = ["id", "type","batchId"]
sim["topology"]["structure"]["data"] = []
for i in range(vertices1.shape[0]):
    sim["topology"]["structure"]["data"].append([i, "A", 0])
for i in range(vertices2.shape[0]):
    sim["topology"]["structure"]["data"].append([i + idOffset2, "A", 1])
for i in range(vertices3.shape[0]):
    sim["topology"]["structure"]["data"].append([i + idOffset3, "A", 1])

sim["topology"]["forceField"] = {}

sim["topology"]["forceField"]["bonds"] = {}
sim["topology"]["forceField"]["bonds"]["type"]       = ["Bond2","Harmonic"]
sim["topology"]["forceField"]["bonds"]["labels"]     = ["id_i","id_j","r0","K"]
sim["topology"]["forceField"]["bonds"]["parameters"] = {}
sim["topology"]["forceField"]["bonds"]["data"]       = []

for bond in bonds1:
    p1 = vertices1[bond[0]]
    p2 = vertices1[bond[1]]
    dst = np.linalg.norm(p1 - p2)
    sim["topology"]["forceField"]["bonds"]["data"].append([int(bond[0]), int(bond[1]), dst, Kb])

for bond in bonds2:
    p1 = vertices2[bond[0]-idOffset2]
    p2 = vertices2[bond[1]-idOffset2]
    dst = np.linalg.norm(p1 - p2)
    sim["topology"]["forceField"]["bonds"]["data"].append([int(bond[0]), int(bond[1]), dst, Kb])

for bond in bonds3:
    p1 = vertices3[bond[0]-idOffset3]
    p2 = vertices3[bond[1]-idOffset3]
    dst = np.linalg.norm(p1 - p2)
    sim["topology"]["forceField"]["bonds"]["data"].append([int(bond[0]), int(bond[1]), dst, Kb])

sim["topology"]["forceField"]["dihedrals"] = {}
sim["topology"]["forceField"]["dihedrals"]["type"]       = ["Bond4","Dihedral"]
sim["topology"]["forceField"]["dihedrals"]["labels"]     = ["id_i","id_j","id_k","id_l","phi0","K","n"]
sim["topology"]["forceField"]["dihedrals"]["parameters"] = {}
sim["topology"]["forceField"]["dihedrals"]["data"]       = []

for d in dihedrals1:
    sim["topology"]["forceField"]["dihedrals"]["data"].append([d[0],d[1],d[2],d[3],0.0,Kd,1])
for d in dihedrals2:
    sim["topology"]["forceField"]["dihedrals"]["data"].append([d[0],d[1],d[2],d[3],0.0,Kd,1])
for d in dihedrals3:
    sim["topology"]["forceField"]["dihedrals"]["data"].append([d[0],d[1],d[2],d[3],0.0,Kd,1])

sim["topology"]["forceField"]["fixed"] = {}
sim["topology"]["forceField"]["fixed"]["type"]       = ["Set1","FixedHarmonicAnisotropicCenterOfMass"]
sim["topology"]["forceField"]["fixed"]["labels"]     = ["idSet_i","K","r0","position"]
sim["topology"]["forceField"]["fixed"]["parameters"] = {}
sim["topology"]["forceField"]["fixed"]["data"]       = [[ids1.copy(), [Kf,Kf,Kf], [0.0,0.0,0.0],[0.0,0.0,L/2]]]

sim["topology"]["forceField"]["comBond"] = {}
sim["topology"]["forceField"]["comBond"]["type"]       = ["Set2","HarmonicBondBetweenCentersOfMass"]
sim["topology"]["forceField"]["comBond"]["labels"]     = ["idSet_i","idSet_j","K","r0"]
sim["topology"]["forceField"]["comBond"]["parameters"] = {}
sim["topology"]["forceField"]["comBond"]["data"]       = [[ids2.copy(), ids3.copy(), Kc, separation]]

#Output

sim["simulationStep"] = {}

sim["simulationStep"]["groups"] = {}
sim["simulationStep"]["groups"]["type"] = ["Groups", "GroupsList"]
sim["simulationStep"]["groups"]["parameters"] = {}
sim["simulationStep"]["groups"]["labels"] = ["name","type","selection"]
sim["simulationStep"]["groups"]["data"] = []
sim["simulationStep"]["groups"]["data"].append(["group1","Ids",ids1.copy()])
sim["simulationStep"]["groups"]["data"].append(["group2","Ids",ids2.copy()])
sim["simulationStep"]["groups"]["data"].append(["group3","Ids",ids3.copy()])

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

sim["simulationStep"]["comPos"] = {}
sim["simulationStep"]["comPos"]["type"] = ["GeometricalMeasure", "CenterOfMassPosition"]
sim["simulationStep"]["comPos"]["parameters"] = {}
sim["simulationStep"]["comPos"]["parameters"]["intervalStep"] = nStepsMeasure
sim["simulationStep"]["comPos"]["parameters"]["outputFilePath"] = "comPos.dat"
sim["simulationStep"]["comPos"]["parameters"]["group"] = "group1"

sim["simulationStep"]["comDst"] = {}
sim["simulationStep"]["comDst"]["type"] = ["GeometricalMeasure", "DistanceBetweenCentersOfMass"]
sim["simulationStep"]["comDst"]["parameters"] = {}
sim["simulationStep"]["comDst"]["parameters"]["intervalStep"] = nStepsMeasure
sim["simulationStep"]["comDst"]["parameters"]["outputFilePath"] = "comDst.dat"
sim["simulationStep"]["comDst"]["labels"] = ["idSet_i","idSet_j"]
sim["simulationStep"]["comDst"]["data"] = [[ids2.copy(), ids3.copy()]]

#Create results folder if it does not exist
if not os.path.exists("results"):
    os.makedirs("results")

sim.write("results/simulation.json")



