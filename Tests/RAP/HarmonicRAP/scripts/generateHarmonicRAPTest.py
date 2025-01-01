import sys,os

import pyUAMMD

import numpy as np
import quaternion as qt
from scipy.spatial.transform import Rotation

import json
import jsbeautifier

from VLMP.utils.geometry.objects.helix import *

#######################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

N  = param["N"]

timeStep = param["timeStep"]

temperature = param["temperature"]

mass   = param["mass"]
radius = param["radius"]

Khrm = param["Khrm"]

Krap = param["Krap"]

helixRadius = param["helixRadius"]
helixPitch  = param["helixPitch"]

viscosity = param["viscosity"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]

perturbation = param["perturbation"]

print("Parameters:")
print(f"N = {N}")
print(f"timeStep = {timeStep}")
print(f"temperature = {temperature}")
print(f"mass = {mass}")
print(f"radius = {radius}")
print(f"Khrm = {Khrm}")
print(f"Krap = {Krap}")
print(f"helixRadius = {helixRadius}")
print(f"helixPitch = {helixPitch}")
print(f"viscosity = {viscosity}")
print(f"nSteps = {nSteps}")
print(f"nStepsOutput = {nStepsOutput}")

#Compute box size
L = 2.0*radius*(N-1)+5.0*2.0*radius
box = [L,L,L]

pos,ori = generateHelix(N, helixRadius, helixPitch, 1.0, 2.0*radius)
# Center the helix
pos -= np.mean(pos, axis=0)

connS,connE = computeConnections(helixRadius, helixPitch, 1.0, 2.0*radius)
Rmatrix     = computeHelixMatrix(helixRadius, helixPitch, 1.0, 2.0*radius)
Rq          = qt.as_float_array(qt.from_rotation_matrix(Rmatrix)).tolist()

#for i in range(N):
#    pos[i] = [i*2.0*radius,0.0,0.0]
#    ori[i] = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
#pos -= np.mean(pos, axis=0)
#
#connS       = np.asarray([ radius,0.0,0.0])
#connE       = np.asarray([-radius,0.0,0.0])
#Rq          = [1.0,0.0,0.0,0.0]
#
#with open("helix.sp","w") as f:
#    # Format is x y z radius type
#    for i in range(N):
#        x,y,z    = pos[i]
#        xS,yS,zS = ori[i]@connS + pos[i]
#        xE,yE,zE = ori[i]@connE + pos[i]
#
#        f.write(f"{x} {y} {z} {radius} 0\n")
#        f.write(f"{xS} {yS} {zS} {radius/4} 1\n")
#        f.write(f"{xE} {yE} {zE} {radius/4} 2\n")

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "RAPTest"

simulation["global"] = {}

simulation["global"]["units"] = {}
simulation["global"]["units"]["type"] = ["Units","None"]

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

simulation["integrator"] = {}

simulation["integrator"]["eulerMaruyamaRigid"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["type"] = ["Brownian", "EulerMaruyamaRigidBody"]
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"] = {}
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["eulerMaruyamaRigid"]["parameters"]["viscosity"] = viscosity

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "eulerMaruyamaRigid", nSteps]
]

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction"]
simulation["state"]["data"] = []

previous = np.zeros(3)
for i in range(N):
    p = pos[i]
    q = qt.from_rotation_matrix(ori[i])

    # Generate small random translation
    translation = np.random.uniform(-perturbation,perturbation,3)
    p += translation

    # Generate small random rotation
    # Generate random rotation axis
    axis = np.random.uniform(-1.0,1.0,3)
    axis /= np.linalg.norm(axis)
    # Generate random rotation angle
    angle = np.random.uniform(-np.pi,np.pi)*perturbation
    rot = Rotation.from_rotvec(angle*axis).as_matrix()
    rot = qt.from_rotation_matrix(rot)

    q = qt.as_float_array(q*rot)
    simulation["state"]["data"].append([i,p.tolist(),q.tolist()])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"]   = []

for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["hrmrap"] = {}
simulation["topology"]["forceField"]["hrmrap"]["type"]   = ["Bond2", "HarmonicRAP"]
simulation["topology"]["forceField"]["hrmrap"]["labels"] = ["id_i", "id_j", "Khrm", "r0", "leftConnection", "rightConnection","Krap", "R"]
simulation["topology"]["forceField"]["hrmrap"]["data"]   = []

for i in range(N-1):
    simulation["topology"]["forceField"]["hrmrap"]["data"].append([i, i+1, Khrm, 0.0, connS.tolist(), connE.tolist(), Krap, Rq])

#Output
simulation["simulationStep"] = {}
simulation["simulationStep"][f"write"] = {}
simulation["simulationStep"][f"write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"][f"write"]["parameters"] = {}
simulation["simulationStep"][f"write"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"][f"write"]["parameters"]["outputFilePath"] = f"output"
simulation["simulationStep"][f"write"]["parameters"]["outputFormat"]   = "itpd"

simulation["simulationStep"][f"writeSPO"] = {}
simulation["simulationStep"][f"writeSPO"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"][f"writeSPO"]["parameters"] = {}
simulation["simulationStep"][f"writeSPO"]["parameters"]["intervalStep"]   = nStepsOutput
simulation["simulationStep"][f"writeSPO"]["parameters"]["outputFilePath"] = f"output"
simulation["simulationStep"][f"writeSPO"]["parameters"]["outputFormat"]   = "spo"

simulation["simulationStep"]["thermo"] = {}
simulation["simulationStep"]["thermo"]["type"] = ["ThermodynamicMeasure","ThermodynamicQuantityMeasure"]
simulation["simulationStep"]["thermo"]["parameters"] = {}
simulation["simulationStep"]["thermo"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["thermo"]["parameters"]["outputFilePath"] = "thermo.dat"

simulation["simulationStep"]["potMeasure"] = {}
simulation["simulationStep"]["potMeasure"]["type"] = ["ParticlesListMeasure","PotentialMeasure"]
simulation["simulationStep"]["potMeasure"]["parameters"] = {
        "intervalStep": 1,
        "outputFilePath": "pot.dat"
        }
simulation["simulationStep"]["potMeasure"]["labels"] = ["id"]
simulation["simulationStep"]["potMeasure"]["data"] = [[i] for i in range(N)]

simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
