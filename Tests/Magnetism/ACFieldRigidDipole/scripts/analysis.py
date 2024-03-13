import sys
from solveFP import computeCycle
import numpy as np
import json
import matplotlib.pyplot as plt

#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

b0            = param["b0"]
f             = param["frequency"]
vis           = param["viscosity"]
kBT           = param["temperature"]
rc            = param["radius"]
stepsPerCycle = param["stepsPerCycle"]
timeStep      = 1/(f*stepsPerCycle)
m0            = 4./3.*np.pi*param["msat"]*rc**3

#Compute the cycles using the FP equation
field_teor, magnetization_teor = computeCycle(kBT, vis, rc, m0, f, b0,
                                              ncycles = 5, pointsPerCycle = 5000,
                                              nindex = 50)

#Read results of the simulation
dataSimulation     = np.loadtxt("./results/output.magnet")
time_simu          = dataSimulation[:,0]
field_simu         = b0*np.sin(2*np.pi*f*time_simu)
magnetization_simu = dataSimulation[:,-1]

#Represent the results
plt.plot(field_teor, magnetization_teor, color = "red", label = "FP")
plt.plot(field_simu, magnetization_simu, marker = "o",
         color = "k", label = "BD", linestyle="None")
plt.xlabel("B (A. U.)")
plt.ylabel("$M_z/M_{max}$")
plt.legend(loc = "upper left")
plt.show()
