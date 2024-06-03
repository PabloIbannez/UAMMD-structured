import json
import JFIO

import numpy as np

def linear(x,y,z):
    return x+y+z

def product(x,y,z):
    return x*y*z

with open("./parameters.json") as f:
    inputData = json.load(f)

tol = 1e-4

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

# Read the simulation data

sim = JFIO.read("./results/simulation.json")

# Read the simulation data
file_path = "./results/potential.dat"
data = np.loadtxt(file_path,skiprows=2)
E_sim = data[:,1]


someError = False
state = sim["state"]["data"]
for particle in state:
    i = int(particle[0])

    x = float(particle[1][0])
    y = float(particle[1][1])
    z = float(particle[1][2])

    p = particle[2]

    E_theo = energyTabulated(x,y,z)*scale*p

    diff = np.abs(E_theo-E_sim[i])
    if diff > tol:
        someError = True
        print("Error for particle ",particle[0]," at position ",x,y,z)
        print("Theoretical energy: ",E_theo," Simulation energy: ",E_sim[0])
        print("Error: ",diff)
        print("")

if not someError:
    print("No errors found in the simulation")



