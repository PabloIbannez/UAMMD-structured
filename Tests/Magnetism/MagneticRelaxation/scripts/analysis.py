import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_function(t, tau):
    return np.exp(-t / tau)

def fit_exponential(data_t, data_y):
    popt, pcov = curve_fit(exp_function, data_t, data_y)
    return popt[0]


#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)


vis = param["viscosity"]
kBT = param["temperature"]
rc = param["radius"]

#-------------------- Rigid dipole approximation test ----------------------#
#Read results of the simulation
dataSimulationRD = np.loadtxt("./results/RelaxationRD.magnet")
timeRD_simu = dataSimulationRD[:,0]
magnetizationRD_simu = dataSimulationRD[:,-1]
t_brown = fit_exponential(timeRD_simu, magnetizationRD_simu)

#-------------------- Neel relaxation ----------------------------#
#Read results of the simulation
dataSimulationF = np.loadtxt("./results/RelaxationFixed.magnet")
timeF_simu = dataSimulationF[:,0]
magnetizationF_simu = dataSimulationF[:,-1]
t_neel = fit_exponential(timeF_simu, magnetizationF_simu)


#-------------------- Combined relaxation ----------------------------#
#Read results of the simulation
dataSimulationC = np.loadtxt("./results/RelaxationCombined.magnet")
timeC_simu = dataSimulationF[:,0]
magnetizationC_simu = dataSimulationC[:,-1]
t_combined = fit_exponential(timeC_simu, magnetizationC_simu)



#------------------------Print relaxation times------------------------------#
t_brown_expected = 4*np.pi*vis*rc**3/kBT
print("Obtained Brown time:", t_brown, "; Expected:", t_brown_expected)
print("Obtained Neel time:", t_neel)
t_combined_expected = 1/(1/t_neel+1/t_brown)
print("Obtained combined relaxation time:", t_combined, "; Expected:", t_combined_expected)

#------------------------Plot relaxation---------------------------------#

magnetizationRD_theory = np.exp(-timeRD_simu/t_brown_expected)
magnetizationC_theory = np.exp(-timeC_simu/t_combined_expected)

#Represent the results
plt.plot(timeRD_simu[0:-1:5], magnetizationRD_simu[0:-1:5], marker = "o",
         color = "red", label = "Simulation Brown", linestyle="None")

plt.plot(timeRD_simu, magnetizationRD_theory,
         color = "red", label = "Theory Brown", linewidth = 2)

plt.plot(timeF_simu[0:-1:5], magnetizationF_simu[0:-1:5], marker = "o",
         color = "green", label = "Simulation Neel", linestyle="None")

plt.plot(timeRD_simu[0:-1:5], magnetizationC_simu[0:-1:5], marker = "o",
         color = "blue", label = "Simulation Neel+Brown", linestyle="None")

plt.plot(timeRD_simu, magnetizationC_theory,
         color = "blue", label = "Theory Neel+Brown", linewidth = 2)

plt.yscale("log")
plt.xlabel("Time (A. U.)")
plt.ylabel("$M_z/M_{max}$")
plt.legend(loc = "lower left")
plt.tight_layout()
plt.show()
