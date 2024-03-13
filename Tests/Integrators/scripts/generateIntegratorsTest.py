import sys,os

import pyUAMMD

import json
import jsbeautifier

with open("parameters.json", "r") as f:
    param = json.load(f)

N       = param["N"]
density = param["density"]

timeStep = param["timeStep"]

temperature = param["temperature"]

mass = param["mass"]
radius = param["radius"]

frictionConstant = param["frictionConstant"]
viscosity = param["viscosity"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]
nStepsMeasure = param["nStepsMeasure"]

#Compute box size
L = (N/density)**(1./3.)
box = [L,L,L]

#Create simulation

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "integratorsTest"

simulation["global"] = {}

simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
simulation["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]

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

simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position", "direction"]
simulation["state"]["data"] = []
for i in range(N):
    simulation["state"]["data"].append([i, [0,0,0], [0,0,0,1.0]])

simulation["topology"] = {}
simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []
for i in range(N):
    simulation["topology"]["structure"]["data"].append([i, "A"])

#Output

simulation["simulationStep"] = {}
simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {}
simulation["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

simulation["simulationStep"]["write"] = {}
simulation["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["write"]["parameters"] = {}
simulation["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
simulation["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
simulation["simulationStep"]["write"]["parameters"]["outputFormat"] = "sp"

simulation["simulationStep"]["MSD_EULER_MARUYAMA"] = {}
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["type"] = ["GeometricalMeasure", "MeanSquareDisplacement"]
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["parameters"] = {}
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["parameters"]["intervalStep"] = nStepsMeasure
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["parameters"]["outputFilePath"] = "msdEulerMaruyama.dat"
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["parameters"]["startStep"] = 0
simulation["simulationStep"]["MSD_EULER_MARUYAMA"]["parameters"]["endStep"]   = nSteps

simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"] = {}
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["type"] = ["GeometricalMeasure", "MeanSquareDisplacement"]
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["parameters"] = {}
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["parameters"]["intervalStep"] = nStepsMeasure
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["parameters"]["outputFilePath"] = "msdEulerMaruyamaRigid.dat"
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["parameters"]["startStep"] = nSteps
simulation["simulationStep"]["MSD_EULER_MARUYAMA_RIGID"]["parameters"]["endStep"]   = 2*nSteps

simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"] = {}
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["type"] = ["GeometricalMeasure", "MeanAngularCorrelation"]
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["parameters"] = {}
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["parameters"]["intervalStep"] = nStepsMeasure
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["parameters"]["outputFilePath"] = "macEulerMaruyamaRigid.dat"
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["parameters"]["startStep"] = nSteps
simulation["simulationStep"]["MAC_EULER_MARUYAMA_RIGID"]["parameters"]["endStep"]   = 2*nSteps

simulation["simulationStep"]["VEL_BBK"] = {}
simulation["simulationStep"]["VEL_BBK"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["VEL_BBK"]["parameters"] = {}
simulation["simulationStep"]["VEL_BBK"]["parameters"]["intervalStep"] = nStepsMeasure
simulation["simulationStep"]["VEL_BBK"]["parameters"]["outputFilePath"] = "BBK"
simulation["simulationStep"]["VEL_BBK"]["parameters"]["outputFormat"]   = "vel"
simulation["simulationStep"]["VEL_BBK"]["parameters"]["startStep"] = 2*nSteps
simulation["simulationStep"]["VEL_BBK"]["parameters"]["endStep"]   = 3*nSteps

simulation["simulationStep"]["VEL_GJF"] = {}
simulation["simulationStep"]["VEL_GJF"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["VEL_GJF"]["parameters"] = {}
simulation["simulationStep"]["VEL_GJF"]["parameters"]["intervalStep"] = nStepsMeasure
simulation["simulationStep"]["VEL_GJF"]["parameters"]["outputFilePath"] = "GJF"
simulation["simulationStep"]["VEL_GJF"]["parameters"]["outputFormat"]   = "vel"
simulation["simulationStep"]["VEL_GJF"]["parameters"]["startStep"] = 2*nSteps
simulation["simulationStep"]["VEL_GJF"]["parameters"]["endStep"]   = 3*nSteps

#Check if ./results folder exists, if not create it
if not os.path.exists("./results"):
    os.makedirs("./results")

simulation.write("./results/simulation.json")
