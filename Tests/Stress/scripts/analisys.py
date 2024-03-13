import numpy as np

import json

def bond_force(p1,p2,K,r0):
    dr = np.asarray(p2) - np.asarray(p1)
    r  = np.linalg.norm(dr)
    return K*(r-r0)*dr/r

def wca_force(p1,p2,eps,sig):
    dr = np.asarray(p2) - np.asarray(p1)
    r  = np.linalg.norm(dr)
    if r > sig:
        return np.zeros(3)
    else:
        return -12*eps*(sig**12/r**13 - sig**6/r**7)*dr/r

def compute_stress(p1,p2,f12):
    dr = np.asarray(p2) - np.asarray(p1)
    s  = -0.5*np.outer(dr,f12)
    #print(s)
    return -0.5*np.outer(dr,f12)

with open("./results/simulation.json", "r") as f:
    sim = json.load(f)

partPositionsData = sim["state"]["data"]

partPositions = []
for p in partPositionsData:
    partPositions.append(p[1])

bondsData = sim["topology"]["forceField"]["bonds"]["data"]

r0 = bondsData[0][2]
K  = bondsData[0][3]

print("K = " , K)
print("r0 = ", r0)

bonds = []
for b in bondsData:
    bonds.append([b[0], b[1]])

vol    = [0.0 for i in range(len(partPositions))]

z = np.asarray([0.0,0.0,0.0])
stress = [compute_stress(z,z,z) for i in range(len(partPositions))]

for b in bonds:
    p1 = partPositions[b[0]]
    p2 = partPositions[b[1]]
    f = bond_force(p1,p2,K,r0)
    stress[b[0]] += compute_stress(p1,p2, f)
    stress[b[1]] += compute_stress(p2,p1,-f)

wcaData = sim["topology"]["forceField"]["wca"]["data"]

eps = wcaData[0][2]
sgm = wcaData[0][3]

print("eps = ", eps)
print("sgm = ", sgm)

for i in range(len(partPositions)):
    for j in range(i+1,len(partPositions)):
        p1 = partPositions[i]
        p2 = partPositions[j]
        f = wca_force(p1,p2,eps,sgm)
        stress[i] += compute_stress(p1,p2, f)
        stress[j] += compute_stress(p2,p1,-f)

for i in range(len(partPositions)):
    invr  = 0.0
    invr2 = 0.0
    for j in range(0,len(partPositions)):
        if i == j:
            continue
        dr = np.asarray(partPositions[i]) - np.asarray(partPositions[j])
        r  = np.linalg.norm(dr)
        invr  += 1.0/r
        invr2 += 1.0/r**2

    radius = 0.5*invr/invr2
    vol[i] = 4.0/3.0*np.pi*radius**3

stressData = []
volData    = []
with open("./results/stress.dat", "r") as f:
    for line in f:
        if "#" not in line:
            x,y,z,r,c,v,xx,xy,xz,yx,yy,yz,zx,zy,zz = line.split()
            st = [xx,xy,xz,yx,yy,yz,zx,zy,zz]
            st = [float(s) for s in st]
            st = np.asarray([st[0:3],st[3:6],st[6:9]])
            stressData.append(st)
            volData.append(float(v))

allOk = True
for i in range(len(partPositions)):
    # Check if zero
    diff = np.isclose(stress[i],stressData[i],atol=1e-5).flatten()
    diff = np.all(diff)
    if not diff:
        print("Stress not equal for particle ", i)
        print("Computed: ", stress[i])
        print("Expected: ", stressData[i])
        allOk = False
    diff = np.isclose(vol[i],volData[i],atol=1e-5)
    if not diff:
        print("Volume not equal for particle ", i)
        print("Computed: ", vol[i])
        print("Expected: ", volData[i])
        allOk = False

if allOk == True:
    print("OK !!!")
else:
    print("Errors found")

