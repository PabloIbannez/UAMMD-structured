import sys
from solveFP import computeCycle
import numpy as np
import json
import matplotlib.pyplot as plt

def computeSelfMobility(L, viscosity, rh):
    a = rh/L;
    a2= a*a;
    a3 = a2*a;
    c = 2.83729
    b = 0.19457;
    a6pref = 16.0*np.pi*np.pi/45.0 + 630.0*b*b;
    return  1.0/(6.0*np.pi*viscosity*rh)*(1.0-c*a+(4.0/3.0)*np.pi*a3-a6pref*a3*a3);

#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

viscosity = param["viscosity"]
timeStep  = param["timeStep"]
lbox      = param["lbox"]
force     = param["force"]

#Read results of the simulation
dataSimulation  = np.loadtxt("./results/pos.sp")
radius          = dataSimulation[0,3]
posz            = dataSimulation[:,2]

time = timeStep*np.arange(len(posz))

selfMobility = computeSelfMobility(lbox, viscosity, radius)
posz_teor = time*selfMobility*force

#Represent the results
plt.plot(time, posz_teor, color = "red", label = "Theory")
plt.plot(time[0:-1:5], posz[0:-1:5], marker = "o", color = "k",
         label = "Simulation", linestyle="None")
plt.xlabel("Time (A. U.)")
plt.ylabel("z (A. U.)")
plt.legend(loc = "upper left")
plt.show()
