import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.integrate import quad

def computeMagnetization(time, gamma0, alpha, b0):
    c = gamma0/(1+alpha*alpha);
    w = c*b0;
    mx = -np.sin(w*time)/np.cosh(w*alpha*time);
    my = np.cos(w*time)/np.cosh(w*alpha*time);
    mz = np.tanh(w*alpha*time);
    return [mx, my, mz];


#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

dt = param["timeStep"]
nsteps = param["nSteps"]
nsave = param["nStepsMeasure"]
b0 = param["b0"]
direction = param["direction"]
field = b0*np.array(direction)
damping = param["damping"]
gyroRatio = param["gyroRatio"]

#Read results of the simulation
dataSimulation = np.loadtxt("./results/output.magnet")

time = np.linspace(0, nsteps*dt, int(nsteps/nsave)+1)
time = time[:-1]

m = computeMagnetization(time, gyroRatio, damping, field[-1])
plt.plot(time, m[0], color = "r", label = "Theory mx")
plt.plot(time, m[1], color = "g", label = "Theory my")
plt.plot(time, m[2], color = "b", label = "Theory mz")
plt.plot(time, dataSimulation[:,0], marker = "o", linestyle="None",
         color = "r", label = "Simulation mx")
plt.plot(time, dataSimulation[:,1], marker = "o", linestyle="None",
         color = "g", label = "Simulation my")
plt.plot(time[0:-1:10], dataSimulation[0:-1:10,2], marker = "o",
         linestyle="None", color = "b", label = "Simulation mz")
plt.xlabel("Time (A. U.)")
plt.ylabel("$M/M_{sat}$")
plt.legend(ncol = 2)
plt.show()
