import os
import json

import numpy as np

import pyUAMMD

def linear(x,y,z):
    return x+y+z

def product(x,y,z):
    return x*y*z

with open("./parameters.json") as f:
    inputData = json.load(f)

sigma  = inputData["sigma"]

Nx = inputData["Nx"]
Ny = inputData["Ny"]
Nz = inputData["Nz"]

N = inputData["N"]

Lx = inputData["Lx"]
Ly = inputData["Ly"]
Lz = inputData["Lz"]

T  = inputData["T"]
dt = inputData["dt"]

nSteps = inputData["nSteps"]

potType = inputData["potType"]
scale   = inputData["scale"]

if potType not in ["linear","product"]:
    print("Potential type not recognized. Exiting.")
    exit()
else:
    if potType == "linear":
        energyTabulated = linear
    else:
        energyTabulated = product

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

diffx = (Lx/Nx)
diffy = (Ly/Ny)
diffz = (Lz/Nz)

x_array = np.linspace(-Lx/2+diffx/2.0,Lx/2-diffx/2.0,Nx)
y_array = np.linspace(-Ly/2+diffy/2.0,Ly/2-diffy/2.0,Ny)
z_array = np.linspace(-Lz/2+diffz/2.0,Lz/2-diffz/2.0,Nz)

X,Y,Z = np.meshgrid(x_array,y_array,z_array,indexing='ij')

E = energyTabulated(X.flatten(),Y.flatten(),Z.flatten())
E = E.reshape(X.shape)
F = np.gradient(E,diffx,diffy,diffz)
F = [-F[0],-F[1],-F[2]]

sim["state"]={}
sim["state"]["labels"]=["id","position"]
sim["state"]["data"]=[]
for i in range(N):
    x,y,z = np.random.rand(3)*np.array([Lx-2*diffx,Ly-2*diffy,Lz-2*diffz])-np.array([(Lx-2*diffx)/2,(Ly-2*diffy)/2,(Lz-2*diffz)/2])
    sim["state"]["data"].append([i,[x,y,z]])

sim["topology"]={}

sim["topology"]["structure"]={}
sim["topology"]["structure"]["labels"]=["id","type"]
sim["topology"]["structure"]["data"]=[]
for i in range(N):
    sim["topology"]["structure"]["data"].append([i,"A"])

sim["topology"]["forceField"]={}

sim["topology"]["forceField"]["ExTab"]={}
sim["topology"]["forceField"]["ExTab"]["type"]=["External","ExternalTabulated"]
sim["topology"]["forceField"]["ExTab"]["parameters"]={"nx":Nx,"ny":Ny,"nz":Nz,"scale":scale}
sim["topology"]["forceField"]["ExTab"]["labels"]=["i","j","k","energy","force"]
sim["topology"]["forceField"]["ExTab"]["data"]=[]
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            sim["topology"]["forceField"]["ExTab"]["data"].append([i,j,k,E[i,j,k],[F[0][i,j,k],F[1][i,j,k],F[2][i,j,k]]])


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
