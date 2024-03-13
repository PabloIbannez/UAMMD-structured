import numpy as np
import sys

#Sums all the pair forces acting on each particle, it should be equal to the total force acting
#on that particle.

threshold     = float(sys.argv[1])
filenameTotal = "results/TotalForces.out"
filenamePair  = "results/PairForces.out"

totalForces = np.loadtxt(filenameTotal)[:,1:]
pairForces  = np.loadtxt(filenamePair)

sumPairForces = np.zeros([4,3])

for line in pairForces:
    sumPairForces[int(line[1]),:] += line[2:]

diff     = np.abs(totalForces-sumPairForces)
maxError = np.max(diff)

if (maxError<threshold):
    print(f"TEST PAIRWISE_FORCES BOND4 PASSED: The maximum error is {maxError}")
else:
    print(f"TEST PAIRWISE_FORCES BOND4 FAILED: The maximum error is {maxError}")
