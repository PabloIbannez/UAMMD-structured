import pyUAMMD
import json5 as json

with open("parameters.json", "r") as f:
    param = json.load(f)

N       = param["N"]
density = param["density"]

timeStep = param["timeStep"]

temperature = param["temperature"]

mass   = param["mass"]
radius = param["radius"]

frictionConstant = param["frictionConstant"]

nSteps        = param["nSteps"]
nStepsOutput  = param["nStepsOutput"]

#Compute box size
L = (N/density)**(1./3.)
box = [L,L,L]

#Create simulation

sim = pyUAMMD.simulation(debug=True)

sim["system"] = {}
sim["system"]["info"] = {}
sim["system"]["info"]["type"] = ["Simulation","Information"]
sim["system"]["info"]["parameters"] = {}
sim["system"]["info"]["parameters"]["name"] = "pythonTest"

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
sim["global"]["types"]["data"]   = [["A", mass, radius, 0.0]]

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
sim["state"]["labels"] = ["id", "position", "direction"]
sim["state"]["data"] = []
for i in range(N):
    sim["state"]["data"].append([i, [0,0,0], [0,0,0,1.0]])

sim["topology"] = {}
sim["topology"]["structure"] = {}
sim["topology"]["structure"]["labels"] = ["id", "type"]
sim["topology"]["structure"]["data"] = []
for i in range(N):
    sim["topology"]["structure"]["data"].append([i, "A"])

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

sim.write("simulation.json",legacy=True)
#sim.run()
