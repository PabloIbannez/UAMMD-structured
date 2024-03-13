import sys
import json

import subprocess

import numpy as np
import matplotlib.pyplot as plt

with open("parameters.json") as f:
    parameters = json.load(f)

nModel = parameters["nModel"]
nSim   = parameters["nSim"]

L    = 25.0
rcut = 10.0
N    = 8000

referenceFile = "scripts/reference/lj_LAMMPS_rdf.dat"

#Compute radial distribution function
rdfs  = {}
index = 0
for i in range(nModel):
    for j in range(nSim):
        fileName = f"results/output{index}.sp"

        print(f"[INFO] Computing radial distribution function for {fileName}")

        #Count number of spanshots
        nSnapshots = 0
        with open(fileName) as f:
            for line in f:
                if line.startswith("#"):
                    nSnapshots += 1

        #Call program "rdf" to compute radial distribution function
        #rdf output0.sp -L 25 -rcut 10 -N 8000 -nbins 1000 -Nsnapshots $nsnap output0.sp > lj_UAMMD_rdf0.dat
        subprocess.call(["rdf", fileName,
                         "-L", str(L), "-rcut", str(rcut),
                         "-N", str(N),
                         "-nbins", "1000",
                         "-Nsnapshots", str(nSnapshots)],
                         stdout=open(f"results/rdf{index}.dat", "w"))

        rdfs[index] = f"results/rdf{index}.dat"

        index += 1

for index,rdf in rdfs.items():
    rdfData = np.loadtxt(rdf)

    r = rdfData[:,0]
    g = rdfData[:,1]

    plt.plot(r,g,'.', label=f"rdf{index}")

#Load reference
rdfRefData = np.loadtxt(referenceFile)

rRef = rdfRefData[:,0]
gRef = rdfRefData[:,1]

plt.plot(rRef,gRef,'-', label="reference")

plt.legend()
plt.show()





