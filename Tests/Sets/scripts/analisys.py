import json

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simpson

kBT = 1.0

def Harmonic(r,K,r0):
    return 0.5*K*(r-r0)*(r-r0)

def boltz_fixed(x,func,funcParams):
    exponent =-func(x,**funcParams)/kBT
    return np.exp(exponent)

def boltz_bond2(x,func,funcParams):
    exponent =-func(x,**funcParams)/kBT
    return x*x*np.exp(exponent)

parametersFile = "parameters.json"

posFile = "./results/comPos.dat"
dstFile = "./results/comDst.dat"

# Read parameters
with open(parametersFile, "r") as f:
    parameters = json.load(f)

Kf = parameters["Kf"] # Fixed point
Kc = parameters["Kc"] # Fixed distance

L = parameters["boxSize"]
R = parameters["radius"]

pos0 = [0.0, 0.0, L/2.0] # Initial position
r0   = 2*R*1.5

# Read data
pos = np.loadtxt(posFile,skiprows=1)
X = pos[:,1]
Y = pos[:,2]
Z = pos[:,3]

fig = plt.figure()

plt.hist(X-pos0[0], bins=50, density=True, label="X",alpha=0.33)
plt.hist(Y-pos0[1], bins=50, density=True, label="Y",alpha=0.33)
plt.hist(Z-pos0[2], bins=50, density=True, label="Z",alpha=0.33)

x = np.linspace(-L/2.0,L/2.0,10000)
prob = boltz_fixed(x,Harmonic,{"K":Kf,"r0":0.0})
Z = simpson(prob,x)

plt.plot(x,prob/Z,label="Harmonic")

plt.xlim(np.min(X-pos0[0]),np.max(X-pos0[0]))

plt.xlabel("Distance")
plt.ylabel("Probability")

# Save figure
plt.savefig("results/posHist.png")
#plt.show()

# Plot distance
dst = np.loadtxt(dstFile,skiprows=1)
dst = dst[:,1]

# Plot histogram of distance
plt.figure()
plt.hist(dst, bins=50, density=True, label="Distance")

x = np.linspace(0.0,L,10000)
prob = boltz_bond2(x,Harmonic,{"K":Kc,"r0":r0})
Z = simpson(prob,x)

plt.plot(x,prob/Z,label="Harmonic")

plt.xlim(np.min(dst),np.max(dst))

plt.xlabel("Distance")
plt.ylabel("Probability")
plt.savefig("results/dstHist.png")
#plt.show()




