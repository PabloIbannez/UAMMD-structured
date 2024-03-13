import numpy as np
import sys, os

fileNumerical  = "results/HessianNumerical.out"
fileAnalytical = "results/HessianAnalytical.out"
folder_name    = os.path.basename(os.getcwd())


threshold = float(sys.argv[1])

dataNumerical  = np.loadtxt(fileNumerical)
dataAnalytical = np.loadtxt(fileAnalytical)

diff = np.abs(dataNumerical-dataAnalytical)
maxError = np.max(diff)/np.max(abs(dataAnalytical))

if maxError > threshold:
    print(f"TEST HESSIAN {folder_name} Failed. The maximum error is {maxError}")
else:
    print(f"TEST HESSIAN {folder_name} Passed. The maximum error is {maxError}")
        
