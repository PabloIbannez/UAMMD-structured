import sys,os

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

from tqdm import tqdm

import json

if len(sys.argv) > 1:
    skip = int(sys.argv[1])
else:
    skip = 10000

def SurfaceLennardJonesType2(z,epsilon,sigma,surfacePos):
    return epsilon*(np.power(sigma/np.abs(z-surfacePos),12)-2.0*np.power(sigma/np.abs(z-surfacePos),6))

def SurfaceWCAType2(z,epsilon,sigma,surfacePos):
    dz = np.abs(z-surfacePos)
    return np.piecewise(dz,[dz<=sigma,dz>sigma],[lambda dz: epsilon*(np.power(sigma/dz,12)-2.0*np.power(sigma/dz,6)) ,lambda dz: 0.0])

def boltz_surf_bound(x,epsilon,sigma,surfacePos,boundEpsilon,boundSigma,boundPosition):
    exponent  = -SurfaceLennardJonesType2(x,epsilon,sigma,surfacePos)/kBT
    exponent += -SurfaceWCAType2(x,boundEpsilon,boundSigma,boundPosition)/kBT
    return np.exp(exponent)

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

boundParams = param["bound"]

boundType      = boundParams["type"]
boundEpsilon   = boundParams["epsilon"]
boundSigma     = boundParams["sigma"]
boundPosition  = boundParams["surfacePosition"]

############################

#List all files in results directory
files = os.listdir("results")

#Loop through files
for file in tqdm(files):
    #If distance in file name
    if "position" in file:
        #Remove extension from file name (. can appear in file names)
        file_name = ".".join(file.split(".")[0:3])

        eps = float(file_name.split("_")[1])
        sgm = float(file_name.split("_")[2])

        data = np.loadtxt("results/"+file,skiprows=skip)[:,2]
        data = data.flatten()

        minDist = 0.0
        maxDist = boundPosition*0.9
        z = np.linspace(minDist,maxDist,1000)

        limInf = sgm*0.5
        limSup = maxDist

        x = np.linspace(limInf,limSup,10000)
        prob = boltz_surf_bound(x,eps,sgm,0.0,boundEpsilon,boundSigma,boundPosition)
        Z = simps(prob,x)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(data, bins=200, density=True,)

        ax.plot(z,boltz_surf_bound(z,eps,sgm,0.0,boundEpsilon,boundSigma,boundPosition)/Z)

        ax.set_xlabel("z")
        ax.set_ylabel("P(z)")
        ax.set_title(f"Probability {eps} {sgm}")
        fig.savefig(f"./results/probability_{eps}_{sgm}.png")
