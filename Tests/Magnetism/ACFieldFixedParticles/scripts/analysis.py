import sys
from solveFP_LLG import computeCycle
import numpy as np
import json
import matplotlib.pyplot as plt

#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

b0        = param["b0"]
f         = param["frequency"]
K         = param["anisotropy"]
damping   = param["damping"]
gyroRatio = param["gyroRatio"]
kBT       = param["temperature"]
rc        = param["radius"]

m0 = 4./3.*np.pi*param["msat"]*rc**3

stepsPerCycle = param["stepsPerCycle"]
timeStep      = 1/(f*stepsPerCycle)

#Compute the cycles using the FP equation
field_teor, magnetization_teor = computeCycle(kBT, rc, K, damping, m0, gyroRatio, f, b0,
                                              ncycles = 10, pointsPerCycle = 5000,
                                              nindex = 30)

#Read results of the simulation
dataSimulation     = np.loadtxt("./results/output.magnet")
time_simu          = dataSimulation[:,0]
field_simu         = b0*np.sin(2*np.pi*f*time_simu)
magnetization_simu = dataSimulation[:,-1]

#Represent the results
plt.plot(field_teor, magnetization_teor, color = "red", label = "FP")
plt.plot(field_simu, magnetization_simu, marker = "o",
         color = "k", label = "LLG", linestyle="None")
plt.xlabel("B (A. U.)")
plt.ylabel("$M_z/M_{max}$")
plt.legend(loc = "upper left")
plt.show()
