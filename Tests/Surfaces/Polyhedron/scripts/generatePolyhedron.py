import os
import json

import numpy as np

import pyUAMMD

with open("./parameters.json") as f:
    inputData = json.load(f)

sigma  = inputData["sigma"]

delta = 1.0*sigma #To avoid WCA to explode

N = inputData["N"]

vertex = inputData["vertex"]
faces  = inputData["faces"]

Lx = inputData["Lx"]
Ly = inputData["Ly"]
Lz = inputData["Lz"]

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

for i in range(int(N)):
    x,y,z = 0*np.random.rand(3)*np.array([Lx,Ly,Lz])-0*np.array([Lx,Ly,Lz])/2
    sim["state"]["data"].append([i,[x,y,z]])

sim["topology"]={}

sim["topology"]["structure"]={}
sim["topology"]["structure"]["labels"]=["id","type"]
sim["topology"]["structure"]["data"]=[]
for i in range(N):
    sim["topology"]["structure"]["data"].append([i,"A"])

sim["topology"]["forceField"]={}

sim["topology"]["forceField"]["cube"]={}
sim["topology"]["forceField"]["cube"]["type"]=["Surface","WCA_Polyhedron"]
sim["topology"]["forceField"]["cube"]["parameters"]={"vertex":vertex, "faces":faces}
sim["topology"]["forceField"]["cube"]["labels"]=["name","epsilon","sigma"]
sim["topology"]["forceField"]["cube"]["data"]=[["A",1,sigma]]

sim["topology"]["forceField"]["p1"]={}
sim["topology"]["forceField"]["p1"]["type"]=["Surface","SurfaceGeneralLennardJonesType2"]
sim["topology"]["forceField"]["p1"]["parameters"]={"surfacePosition":-Lz/2}
sim["topology"]["forceField"]["p1"]["labels"]=["name","epsilon","sigma"]
sim["topology"]["forceField"]["p1"]["data"]=[["A",1,sigma]]

sim["topology"]["forceField"]["p2"]={}
sim["topology"]["forceField"]["p2"]["type"]=["Surface","SurfaceGeneralLennardJonesType2"]
sim["topology"]["forceField"]["p2"]["parameters"]={"surfacePosition":Lz/2}
sim["topology"]["forceField"]["p2"]["labels"]=["name","epsilon","sigma"]
sim["topology"]["forceField"]["p2"]["data"]=[["A",1,sigma]]

sim["simulationStep"]={}

sim["simulationStep"]["info"]={}
sim["simulationStep"]["info"]["type"]=["ParticlesListMeasure","PotentialMeasure"]
sim["simulationStep"]["info"]["parameters"]={
        "intervalStep":nSteps,
        "outputFilePath":"potential.dat"
        }
sim["simulationStep"]["info"]["labels"]=["id"]
sim["simulationStep"]["info"]["data"]=[[i] for i in range(N)]


#sim["simulationStep"]["saveState"]={}
#sim["simulationStep"]["saveState"]["type"]=["WriteStep","WriteStep"]
#sim["simulationStep"]["saveState"]["parameters"]={
#    "intervalStep":np.ceil(nSteps/100),
#    "outputFilePath":"output",
#    "outputFormat":"sp"
#    }

try:
    os.makedirs("results")
except:
    pass

sim.write("results/simulation.json")
