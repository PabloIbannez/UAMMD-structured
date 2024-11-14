import pyUAMMD
import numpy as np


def distributeRandomlyParticles(L, nParticles):
    particles = []
    pos = np.zeros((nParticles, 3))
    for i in range(nParticles):
        x = np.random.uniform(-0.5*L, 0.5*L)
        y = np.random.uniform(-0.5*L, 0.5*L)
        z = np.random.uniform(-0.5*L, 0.5*L)
        pos[i,:] = np.array([x,y,z])
    return pos

nParticles        = 1000
L                 = 100
radius            = 1.0

numberOfSteps     = 100000
numberStepsOutput = 1000
timeStep          = 0.01

temperature       = 1.0
viscosity         = 1.0
psi               = 0.6
tolerance         = 1e-3

box = [L,L,L]
pos = distributeRandomlyParticles(L, nParticles)

simulation = pyUAMMD.simulation()

simulation["system"] = {}
simulation["system"]["info"] = {}
simulation["system"]["info"]["type"] = ["Simulation","Information"]
simulation["system"]["info"]["parameters"] = {}
simulation["system"]["info"]["parameters"]["name"] = "ExamplePSE"

simulation["global"] = {}
simulation["global"]["types"] = {}
simulation["global"]["types"]["type"]   = ["Types","Basic"]
simulation["global"]["types"]["labels"] = ["name", "radius", "mass", "charge"]
simulation["global"]["types"]["data"]  = [["A", radius, 0, 0]]

simulation["global"]["ensemble"] = {}
simulation["global"]["ensemble"]["type"]   = ["Ensemble","NVT"]
simulation["global"]["ensemble"]["labels"] = ["box", "temperature"]
simulation["global"]["ensemble"]["data"]   = [[box, temperature]]


simulation["integrator"] = {}
simulation["integrator"]["PSE"] = {}
simulation["integrator"]["PSE"]["type"] = ["BDHITriplePeriodic", "PositivelySplitEwald"]
simulation["integrator"]["PSE"]["parameters"] = {}
simulation["integrator"]["PSE"]["parameters"]["timeStep"] = timeStep
simulation["integrator"]["PSE"]["parameters"]["viscosity"] = viscosity
#simulation["integrator"]["PSE"]["parameters"]["temperature"] = temperature
simulation["integrator"]["PSE"]["parameters"]["hydrodynamicRadius"] = radius
simulation["integrator"]["PSE"]["parameters"]["tolerance"] = tolerance
simulation["integrator"]["PSE"]["parameters"]["psi"] = psi

simulation["integrator"]["schedule"] = {}
simulation["integrator"]["schedule"]["type"] = ["Schedule", "Integrator"]
simulation["integrator"]["schedule"]["labels"] = ["order", "integrator","steps"]
simulation["integrator"]["schedule"]["data"] = [
    [1, "PSE", numberOfSteps],
]


simulation["state"] = {}
simulation["state"]["labels"] = ["id", "position"]
simulation["state"]["data"] = []

for i,p in enumerate(pos):
    simulation["state"]["data"] += [[i, p.tolist()]]

simulation["topology"] = {}

simulation["topology"]["structure"] = {}
simulation["topology"]["structure"]["labels"] = ["id", "type"]
simulation["topology"]["structure"]["data"] = []

simulation["topology"]["forceField"] = {}
simulation["topology"]["forceField"]["forceCte"]   = {}
simulation["topology"]["forceField"]["forceCte"]["type"] = ["External","ConstantForce"]
simulation["topology"]["forceField"]["forceCte"]["parameters"] = {"constantForce":[0,0,1]}

for i in range(len(pos)):
    simulation["topology"]["structure"]["data"] += [[i, "A"]]


simulation["simulationStep"] = {}
simulation["simulationStep"]["info"] = {}
simulation["simulationStep"]["info"]["type"] = ["UtilsStep", "InfoStep"]
simulation["simulationStep"]["info"]["parameters"] = {"intervalStep":numberStepsOutput}

simulation["simulationStep"]["pos"] = {}
simulation["simulationStep"]["pos"]["type"] = ["WriteStep", "WriteStep"]
simulation["simulationStep"]["pos"]["parameters"] = {"intervalStep":numberStepsOutput,
                                                     "outputFilePath": "output",
                                                     "outputFormat":"sp"}
simulation.write(f"simulation.json")
