import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.integrate import quad

def probability_ctheta(ctheta, K, volume, temperature, nPoints):
    KV_div_temp = K*volume/temperature
    return np.exp(-KV_div_temp)*np.exp(KV_div_temp*ctheta**2)

def integrate_probability_ctheta(K, volume, temperature, nPoints):
    integral, _ = quad(probability_ctheta, -1, 1, args=(K, volume, temperature, nPoints))
    return integral

def computeProbabilityDistribution(K, volume, temperature, nPoints = 1000):
    normalization = integrate_probability_ctheta(K, volume, temperature, nPoints)
    ctheta = np.linspace(-1,1,nPoints)
    prob = probability_ctheta(ctheta, K, volume, temperature, nPoints)
    return ctheta, prob/normalization

#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

temperature = param["temperature"]
volume = param["volume"]
anisotropy = param["anisotropy"]

#Read results of the simulation
dataSimulation = np.loadtxt("./results/output.magnet")

#Compute the cosine of the direction of the magnetization with the easy axes
costheta = dataSimulation[:,-1]/np.sqrt(np.dot(dataSimulation[0,:],dataSimulation[0,:]))

#Compute the Boltzman distribution of cos(\theta) using E = KV*sin^2(\theta).
#p(cos(\theta)) = exp(-E/kBT)/Z
ctheta_B, p_B = computeProbabilityDistribution(anisotropy, volume, temperature)

#Plot the results
plt.hist(costheta, bins = 50, density = True, label = "Simulation")
plt.plot(ctheta_B, p_B, label = "Boltzman-Distribution")
plt.xlabel(r"$cos(\theta)$")
plt.ylabel(r"$p(cos(\theta))$")
plt.show()
