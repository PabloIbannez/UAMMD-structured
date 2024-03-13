import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import json
from scipy.integrate import nquad

#Not exact!!!
def computeTheoreticalAngle(permeability, m0, viscosity, hydroRadius, r0, time):
    t1 = -permeability*m0**2*time/(16*np.pi**2*viscosity*hydroRadius**3*r0**3)
    return 2*np.arctan(np.exp(t1))

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
rij          = param["rij"][-1]
hydroRadius  = (3*volume/(4*np.pi))**(1./3.)
#Read results of the simulation
dataSimulation = np.loadtxt("./results/output.magnet")
m0 = msat*volume
ctheta = dataSimulation[1::2, 2]/m0
time = np.arange(len(ctheta))*timeStep*nsave
thetaTheo = computeTheoreticalAngle(permeability, m0, viscosity, hydroRadius, rij, time)
plt.plot(time[:-1:8], ctheta[:-1:8], label = "Simulation",
         marker = "o", color = "k", linestyle="None")
plt.plot(time, np.cos(thetaTheo), label = "Theory", color = "red")
plt.xlabel("Time (A. U.)")
plt.ylabel(r"$cos(\theta)$")
plt.legend()
plt.show()
