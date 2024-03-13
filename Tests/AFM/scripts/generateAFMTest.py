import sys,os
import json5 as json

import pyUAMMD

with open("parameters.json", "r") as f:
    param = json.load(f)

nAFM    = param["nAFM"]

distanceFactor = param["distanceFactor"]

N       = param["N"]
tipPos  = param["tipPos"]
chipPos = param["chipPos"]
epsilon = param["epsilon"]
sigma   = param["sigma"]
K       = param["K"]
Kxy     = param["Kxy"]
tipVel  = param["tipVelocity"]

box     = param["box"]

timeStep = param["timeStep"]

temperature = param["temperature"]

massTip   = param["massTip"]
radiusTip = param["radiusTip"]

massSample   = param["massSample"]
radiusSample = param["radiusSample"]

frictionConstant = param["frictionConstant"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]

#Create simulation

simMerged = None

for i in range(nAFM):

    sim = pyUAMMD.simulation(debug=True)

    sim["system"] = {}
    sim["system"]["info"] = {}
    sim["system"]["info"]["type"] = ["Simulation","Information"]
    sim["system"]["info"]["parameters"] = {}
    sim["system"]["info"]["parameters"]["name"] = "AFMTest"

    sim["global"] = {}

    sim["global"]["types"] = {}
    sim["global"]["types"]["type"]   = ["Types","Basic"]
    sim["global"]["types"]["labels"] = ["name", "mass", "radius", "charge"]
    sim["global"]["types"]["data"]   = [
        ["TIP", massTip, radiusTip, 0.0],
        ["SAMPLE", massSample, radiusSample, 0.0]
    ]

    sim["global"]["ensemble"] = {}
    sim["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
    sim["global"]["ensemble"]["labels"] = ["box", "temperature"]
    sim["global"]["ensemble"]["data"]   = [[box, temperature]]

    sim["integrator"] = {}

    sim["integrator"]["bbk"] = {}
    sim["integrator"]["bbk"]["type"] = ["Langevin", "BBK"]
    sim["integrator"]["bbk"]["parameters"] = {}
    sim["integrator"]["bbk"]["parameters"]["timeStep"] = timeStep
    sim["integrator"]["bbk"]["parameters"]["frictionConstant"] = frictionConstant

    sim["integrator"]["schedule"] = {}
    sim["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
    sim["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
    sim["integrator"]["schedule"]["data"] = [
        [1, "bbk", nSteps]
    ]

    sim["state"] = {}
    sim["state"]["labels"] = ["id", "position"]
    sim["state"]["data"] = []
    sim["state"]["data"].append([0, [0,0,tipPos[i]]])
    for j in range(N[i]):
        sim["state"]["data"].append([j+1, [0,0,tipPos[i]-distanceFactor*(radiusTip+radiusSample)]])

    sim["topology"] = {}
    sim["topology"]["structure"] = {}
    sim["topology"]["structure"]["labels"] = ["id", "type"]
    sim["topology"]["structure"]["data"] = []
    sim["topology"]["structure"]["data"].append([0, "TIP"])
    for j in range(N[i]):
        sim["topology"]["structure"]["data"].append([j+1, "SAMPLE"])

    #Output

    sim["topology"]["forceField"] = {}
    sim["topology"]["forceField"]["AFM"] = {}
    sim["topology"]["forceField"]["AFM"]["parameters"] = {}

    sim["topology"]["forceField"]["AFM"]["type"] = ["AFM", "SphericalTip"]
    sim["topology"]["forceField"]["AFM"]["labels"] = ["idSet_i", "idSet_j", "epsilon","sigma","K","Kxy","tipVelocity","startChipPosition","indentationStartStep"]
    sim["topology"]["forceField"]["AFM"]["data"] = []
    sim["topology"]["forceField"]["AFM"]["data"].append([[0],list(range(1,N[i]+1)),
                                                         epsilon[i],
                                                         sigma[i],
                                                         K[i],Kxy[i],
                                                         tipVel[i],
                                                         [0.0,0.0,chipPos[i]],
                                                         0])

    sim["simulationStep"] = {}
    sim["simulationStep"]["info"] = {}
    sim["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
    sim["simulationStep"]["info"]["parameters"] = {}
    sim["simulationStep"]["info"]["parameters"]["intervalStep"] = nStepsOutput

    sim["simulationStep"]["tip"] = {}
    sim["simulationStep"]["tip"]["type"] = ["ParticlesListMeasure", "PotentialMeasure"]
    sim["simulationStep"]["tip"]["parameters"] = {}
    sim["simulationStep"]["tip"]["parameters"]["intervalStep"] = nStepsOutput
    sim["simulationStep"]["tip"]["parameters"]["outputFilePath"] = "tip.dat"
    sim["simulationStep"]["tip"]["labels"] = ["id"]
    sim["simulationStep"]["tip"]["data"]   = [[0]]

    sim["simulationStep"]["write"] = {}
    sim["simulationStep"]["write"]["type"] = ["WriteStep", "WriteStep"]
    sim["simulationStep"]["write"]["parameters"] = {}
    sim["simulationStep"]["write"]["parameters"]["intervalStep"] = nStepsOutput
    sim["simulationStep"]["write"]["parameters"]["outputFilePath"] = "output"
    sim["simulationStep"]["write"]["parameters"]["outputFormat"] = "sp"

    if simMerged == None:
        simMerged = sim
    else:
        simMerged.append(sim)


# Create results folder if it does not exist yet
if not os.path.exists('results'):
    os.mkdir('results')

simMerged.write("results/simulation.json")
