import sys,os

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

import json

if len(sys.argv) > 1:
    skip = int(sys.argv[1])
else:
    skip = 300

def LennardJonesType2(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: epsilon*(np.power(sigma/r,12)-2.0*np.power(sigma/r,6)),lambda r: 0.0])

def Harmonic(r,K,r0):
    return 0.5*K*(r-r0)*(r-r0)

def MaxDistanceRestraint(r,K,maxDistance):
    return np.piecewise(r,[r<=maxDistance,r>maxDistance],[lambda r: 0.0,lambda r: Harmonic(r,K,maxDistance)])

def boltz_bond2_bound(x,func,funcParams,boundK,boundMaxDistance):
    exponent =-func(x,**funcParams)/kBT
    exponent += -MaxDistanceRestraint(x,boundK,boundMaxDistance)/kBT
    return x*x*np.exp(exponent)

############################

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

units = param["units"]
if units == "none":
    kB      = 1.0
    ELECOEF = 1.0
elif units == "KcalMol_A":
    kB = 1.987191E-03
    ELECOEF = 332.0716

kBT = kB*param["temperature"]

boundParams = param["bond2bound"]

boundType        = boundParams["type"]
boundK           = boundParams["K"]
boundMaxDistance = boundParams["distance"]

############################

#List all files in results directory
files = os.listdir("results")

#Loop through files
for file in files:
    #If distance in file name
    if "distance" in file:
        #Remove extension from file name (. can appear in file names)
        file_name = ".".join(file.split(".")[0:3])

        eps = float(file_name.split("_")[1])
        sgm = float(file_name.split("_")[2])

        data = np.loadtxt("results/"+file,skiprows=skip)
        data = data.flatten()

        minDist = 0.0
        maxDist = boundMaxDistance+2.5
        r = np.linspace(minDist,maxDist,1000)

        limInf = sgm*0.5
        limSup = boundMaxDistance*2.0

        x = np.linspace(limInf,limSup,10000)
        prob = boltz_bond2_bound(x,LennardJonesType2,{"epsilon":eps,"sigma":sgm,"cutOffFactor":2.5},boundK,boundMaxDistance)
        Z = simps(prob,x)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(data, bins=200, density=True,)

        ax.plot(r,boltz_bond2_bound(r,LennardJonesType2,{"epsilon":eps,"sigma":sgm,"cutOffFactor":2.5},boundK,boundMaxDistance)/Z,)

        ax.set_xlabel("r")
        ax.set_ylabel("P(r)")
        ax.set_title(f"Probability {eps} {sgm}")
        fig.savefig(f"./results/probability_{eps}_{sgm}.png")
