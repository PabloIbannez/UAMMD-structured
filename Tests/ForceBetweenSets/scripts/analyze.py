import numpy as np
import sys

filenameSets  = "results/SetsForces.out"
filenamePair  = "results/PairForces.out"

setsForces  = np.loadtxt(filenameSets)
pairForces  = np.loadtxt(filenamePair)

f_oddEven = 0
f_evenOdd = 0

threshold = float(sys.argv[1])
for i, force in enumerate(pairForces):
    if force[0] % 2 == 0:
        f_oddEven += force[2:]
    else:
        f_evenOdd += force[2:]

error_oddEven = abs((f_oddEven-setsForces[0,2:])/setsForces[0,2:])
error_evenOdd = abs((f_evenOdd-setsForces[1,2:])/setsForces[1,2:])

maxError = max(max(error_oddEven), max(error_evenOdd))

if maxError > threshold:
    print(f"TEST FAILED: The maximum error is {maxError}")
else:
    print(f"TEST PASSED: The maximum error is {maxError}")
    
