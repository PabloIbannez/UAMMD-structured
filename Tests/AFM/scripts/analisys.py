import os
import json
import numpy as np

with open("parameters.json", "r") as f:
    parameters = json.load(f)

def harmonic_force_Z(dz, z0, k):
    return k * (z0 - dz)

def generalTip_force_Z(dz, sigma, epsilon,tipRadius):

    if epsilon < 0:
        dst = abs(dz)
        eps = abs(epsilon)

        A =  1.0
        B = -2.0
    else:
        dst = abs(dz)

        A =  0.0
        B =  1.0

    r2 = (dst - tipRadius)**2

    invr2   = 1.0 / r2
    sinvr2  = (sigma * sigma) * invr2
    sinvr6  = sinvr2 ** 3
    sinvr12 = sinvr6 ** 2

    fmod = eps * 6.0 * ( 2.0 * A * sinvr12 + B * sinvr6 )/abs(dst-tipRadius)

    return fmod

nAFM = parameters["nAFM"]

distanceFactor = parameters["distanceFactor"]

N       = parameters["N"]
tipPos  = parameters["tipPos"]
chipPos = parameters["chipPos"]
epsilon = parameters["epsilon"]
sigma   = parameters["sigma"]
K       = parameters["K"]

radiusTip    = parameters["radiusTip"]
radiusSample = parameters["radiusSample"]

resultsFile = "results/tip.dat"
#Check if the file exists
if not os.path.isfile(resultsFile):
    print("File not found: ", resultsFile)
    exit()

results = np.loadtxt(resultsFile)
resultsDict = {}
for r in results:
    resultsDict[int(r[0])] = r[4]

for i in range(nAFM):

    if i > 0:
        tipId = N[i-1] + 1 + tipId
    else:
        tipId = 0

    dz = tipPos[i] - chipPos[i]
    k  = K[i]

    forceAFM = harmonic_force_Z(dz, 0, k)
    forceParticles = 0.0
    for j in range(N[i]):
        pz = np.abs(tipPos[i] - distanceFactor*(radiusTip+radiusSample))
        dz = pz - tipPos[i] #Force over the tip
        forceParticles += generalTip_force_Z(dz, sigma[i], epsilon[i], radiusTip)

    forceTotal = forceAFM + forceParticles

    print(f"Expected force for AFM {i}: {forceAFM:.2f} (AFM), {forceParticles:.2f} (Particles), {forceTotal} (Total). Simulation: {resultsDict[tipId]:.2f}")



