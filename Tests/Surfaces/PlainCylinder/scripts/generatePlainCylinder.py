import os
import json

import numpy as np

import pyUAMMD

with open("./parameters.json") as f:
    inputData = json.load(f)

sigma  = inputData["sigma"]

delta = 1.0*sigma #To avoid WCA to explode

N = inputData["N"]

Lx = inputData["Lx"]
Ly = inputData["Ly"]
Lz = inputData["Lz"]

plainPosition = inputData["plainPosition"]
cylinderRadius = inputData["cylinderRadius"]

T  = inputData["T"]
dt = inputData["dt"]

nSteps = inputData["nSteps"]

sim=pyUAMMD.simulation()

sim["system"]={}
sim["system"]["info"] = {}
sim["system"]["info"]["type"] = ["Simulation","Information"]
sim["system"]["info"]["parameters"] = {"name":"IdealGas_PlainCyl"}


sim["global"]={}

sim["global"]["fundamental"] = {"type":["Fundamental","Time"]}

sim["global"]["ensemble"]={}
sim["global"]["ensemble"]["labels"] = ["box","temperature"]
sim["global"]["ensemble"]["data"] = [[[Lx,Ly,Lz],T]]
sim["global"]["ensemble"]["type"] = ["Ensemble","NVT"]

sim["global"]["types"] = {}
sim["global"]["types"]["labels"] =  ["name", "mass", "radius", "charge"]
sim["global"]["types"]["data"] = [["A",1,sigma,0.0]]
sim["global"]["types"]["type"] = ["Types","Basic"]

sim["integrator"]={}

sim["integrator"]["schedule"]={}
sim["integrator"]["schedule"]["type"]=["Schedule","Integrator"]
sim["integrator"]["schedule"]["data"]=[[1,"verlet",nSteps]]
sim["integrator"]["schedule"]["labels"]=["order","integrator","steps"]

sim["integrator"]["verlet"]={}
sim["integrator"]["verlet"]["parameters"]={}
sim["integrator"]["verlet"]["parameters"]["timeStep"] = dt
sim["integrator"]["verlet"]["type"]=["Verlet","VelocityVerlet"]

sim["state"]={}
sim["state"]["labels"]=["id","position"]
sim["state"]["data"]=[]

for i in range(int(N/2)):
    x,y,z = np.random.rand(3)*np.array([Lx-delta,Ly-delta,Lz-delta])+delta
    sim["state"]["data"].append([i,[x,y,z]])

for i in range(int(N/2),N):
    r = np.random.rand()*(cylinderRadius-delta)
    theta = np.random.rand()*2*np.pi
    z = -np.random.rand()*Lz
    x = r*np.sin(theta)
    y = r*np.cos(theta)
    sim["state"]["data"].append([i,[x,y,z]])

sim["topology"]={}

sim["topology"]["structure"]={}
sim["topology"]["structure"]["labels"]=["id","type"]
sim["topology"]["structure"]["data"]=[]
for i in range(N):
    sim["topology"]["structure"]["data"].append([i,"A"])

sim["topology"]["forceField"]={}

sim["topology"]["forceField"]["wca"]={}
sim["topology"]["forceField"]["wca"]["type"]=["Surface","WCA_PlainCylinder"]
sim["topology"]["forceField"]["wca"]["parameters"]={"plainPosition":plainPosition,"cylinderRadius":cylinderRadius}
sim["topology"]["forceField"]["wca"]["labels"]=["name","epsilon","sigma"]
sim["topology"]["forceField"]["wca"]["data"]=[["A",1,sigma]]

sim["simulationStep"]={}

sim["simulationStep"]["info"]={}
sim["simulationStep"]["info"]["type"]=["ParticlesListMeasure","PotentialMeasure"]
sim["simulationStep"]["info"]["parameters"]={
        "intervalStep":1,
        "outputFilePath":"potential.dat"
        }
sim["simulationStep"]["info"]["labels"]=["id"]
sim["simulationStep"]["info"]["data"]=[[i] for i in range(N)]


sim["simulationStep"]["saveState"]={}
sim["simulationStep"]["saveState"]["type"]=["WriteStep","WriteStep"]
sim["simulationStep"]["saveState"]["parameters"]={
    "intervalStep":1,
    "outputFilePath":"output",
    "outputFormat":"xyz"
    }

try:
    os.makedirs("results")
except:
    pass

sim.write("results/simulation.json")
