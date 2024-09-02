import json
import JFIO

import numpy as np

with open("./parameters.json") as f:
    inputData = json.load(f)

tol = 1e-7

sigma  = inputData["sigma"]
epsilon = inputData["epsilon"]

z0 = inputData["plainPosition"]
Rcyl = inputData["cylinderRadius"]

N = inputData["N"]

Lx = inputData["Lx"]
Ly = inputData["Ly"]
Lz = inputData["Lz"]

T  = inputData["T"]
dt = inputData["dt"]

nSteps = inputData["nSteps"]

# Read the simulation data

sim = JFIO.read("./results/simulation.json")

# Read the simulation data
file_path = "./results/potential.dat"
data = np.loadtxt(file_path,skiprows=2)
E_sim = data[:,1]

def gradientDistance(x,y,z):
                                
    r =  np.sqrt( x* x+ y* y); #//This Surface is defined in Cylindric coordinates

    r_cyl   = Rcyl - r;
    r_plain =  z-z0;

    if ( z<z0 and r<Rcyl): #//Under the plain Inside the cylinder
        r_wall_square = r_cyl*r_cyl;

    elif ( z>z0 and r>Rcyl): #//Above the plain outside the cylinder
        r_wall_square = r_plain*r_plain;

    elif ( z>z0 and r<Rcyl): #//Above the plain insider the cylinder (corner)
        r_wall_square = (r_cyl*r_cyl+r_plain*r_plain);
        rwall = np.sqrt(r_wall_square);

    else:
        return 0.0

    return r_wall_square

def energy(x,y,z):

    gradDistance = gradientDistance(x,y,z);
    r2wall = 0
    r2wall = gradDistance*1

    if(r2wall > sigma*sigma*1.259921):
        return 0.0;

    if(r2wall == 0.0):
        return np.nan

    invr2wall = sigma*sigma/r2wall;
    invr6wall = invr2wall*invr2wall*invr2wall;

    e = 4.0*epsilon*(invr6wall*invr6wall-invr6wall)+epsilon;

    return e



someError = False
state = sim["state"]["data"]
for particle in state:
    i = int(particle[0])

    x = float(particle[1][0])
    y = float(particle[1][1])
    z = float(particle[1][2])

    E_theo = energy(x,y,z)

    diff = np.abs(E_theo-E_sim[i])
    if diff > tol:
        someError = True
        print("Error for particle ",particle[0]," at  tion ",x,y,z)
        print("Theoretical energy: ",E_theo," Simulation energy: ",E_sim[i])
        print("Error: ",diff)
        print("")

if not someError:
    print("No errors found in the simulation")



