import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import json
from scipy.integrate import nquad


def computeTheoreticalDistance(permeability, m0, viscosity, hydroRadius, r0, time):
    return (-10*permeability*m0**2/(4*np.pi**2*viscosity*hydroRadius)*time+r0**5)**(1./5.)

#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

temperature  = param["temperature"]
volume       = param["volume"]
msat         = param["msat"]
permeability = param["permeability"]
timeStep     = param["timeStep"]
nsave        = param["nStepsOutput"]
viscosity    = param["viscosity"]
hydroRadius  = (3*volume/(4*np.pi))**(1./3.)
#Read results of the simulation
dataSimulation = np.loadtxt("./results/output.spm")
m0 = msat*volume
posi = dataSimulation[::4, 2]
posj = dataSimulation[2::4, 2]
rij = posj-posi
time = np.arange(len(rij))*timeStep*nsave
rijTheo = computeTheoreticalDistance(permeability, m0, viscosity, hydroRadius, rij[0], time)
plt.plot(time[:-1:8], rij[:-1:8], label = "Simulation",
         marker = "o", color = "k", linestyle="None")
plt.plot(time, rijTheo, label = "Theory", color = "red")
plt.xlabel("Time (A. U.)")
plt.ylabel("r (A. U.)")
plt.legend()
plt.show()
