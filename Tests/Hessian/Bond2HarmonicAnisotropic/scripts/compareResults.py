import numpy as np
import matplotlib.pyplot as plt
import json


#Read data used in the simulation
with open("parameters.json", "r") as f:
    param = json.load(f)

radius    = param["radius"]
vwallL    = param["vwall"]
vwall     = complex(vwallL[0], vwallL[1])
k         = param["k"]

iterationsList = param["iterationsList"]

with open("results/ImpedanceFromMobility.out", "r") as f:
    resultsMobility = json.load(f)

impedanceRealMob = resultsMobility["graphics"]["y1"][0]
impedanceImagMob = resultsMobility["graphics"]["y2"][0]

impedanceRealIter = []
impedanceImagIter = []
impedanceRealIterb = []
impedanceImagIterb = []

for nIter in iterationsList:
    filename = "results/Impedanceb"+str(nIter)+".out"
    with open(filename, "r") as f:
        resultsMobility = json.load(f)
        impedanceRealIterb += resultsMobility["graphics"]["y0"]
        impedanceImagIterb += resultsMobility["graphics"]["y1"]



Lxy    = param["Lxy"]
theta = np.pi*radius**2/Lxy**2
plt.plot(iterationsList, np.array(impedanceRealIterb)/theta, ".-", color = "red", label = "$Re(Z_{iter})$")
plt.plot(iterationsList, np.array(impedanceImagIterb)/theta, ".-", color = "blue", label = "$Im(Z_{iter})$")

plt.plot(iterationsList, len(iterationsList)*[impedanceRealMob/theta], color = "red", label = "$Re(Z_{mob})$")
plt.plot(iterationsList, len(iterationsList)*[impedanceImagMob/theta], color = "blue", label = "$Im(Z_{mob})$")

print("Z real: ", impedanceRealMob/theta)
print("Z imag: ", impedanceImagMob/theta)


plt.xlabel("Iteration")
plt.ylabel("$Z/Z_{fluid}$")
plt.legend(fontsize = 20, loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("results/compareImpedance.pdf", bbox_inches='tight')
plt.close()

errorRe = abs(np.array(impedanceRealIterb)-impedanceRealMob)/abs(impedanceRealMob)
errorIm = abs(np.array(impedanceImagIterb)-impedanceImagMob)/abs(impedanceImagMob)

plt.plot(iterationsList, errorRe, ".", color = "red", label = "$Re$")
plt.plot(iterationsList, errorIm, ".", color = "blue", label = "$Im$")
plt.yscale("log")
plt.xlabel("Iteration")
plt.ylabel("Relative error")
plt.show()
