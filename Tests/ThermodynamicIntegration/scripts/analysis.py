import numpy as np
import json
import sympy as sp

with open('parameters.json') as f:
    param = json.load(f)

lambdaVal  = param["lambda"]
epsilonVal = param["epsilon"]
sigmaVal   = param["sigma"]
alphaVal   = param["alpha"]
nVal       = param["n"]

###########################################


# Define the symbolic variables for the coordinates of the two points
x1, y1, z1, x2, y2, z2, r = sp.symbols('x1 y1 z1 x2 y2 z2 r')

# Define the distance between the two points
dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
r = sp.sqrt(dx**2 + dy**2 + dz**2)

K, r0, epsilon, sigma, lamb, alpha = sp.symbols('K r0 epsilon sigma lambda alpha')

lambHarmonic = 0.5*K*lamb**nVal*(r - r0)**2
lambFixed    = 0.5*K*lamb**nVal*(dx-r0)**2+0.5*K*lamb**nVal*(dy-r0)**2+0.5*K*lamb**nVal*(dz-r0)**2

den = alpha*(1.0-lamb)**2+(r/sigma)**6
attractive = lamb**nVal/den
repulsive  = lamb**nVal/den**2

lj = epsilon*(repulsive-2.0*attractive)

def bondsEnergyCompute(ids,pos,r0val,Kval,lambdaVal):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    return lambHarmonic.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

def bondsForceCompute(cid,ids,pos,r0val,Kval,lambdaVal):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    dxVal = x2val - x1val
    dyVal = y2val - y1val
    dzVal = z2val - z1val

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    fx = sp.diff(lambHarmonic,x1).subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val),(r0,r0val),(K,Kval),(lamb,lambdaVal)])
    fy = sp.diff(lambHarmonic,y1).subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val),(r0,r0val),(K,Kval),(lamb,lambdaVal)])
    fz = sp.diff(lambHarmonic,z1).subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

    f = np.array([float(fx),float(fy),float(fz)])
    if cid == id2:
        return -f
    return f

def bondsLambdaDerivativeCompute(cid,ids,pos,r0val,Kval,lambdaVal):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    dxVal = x2val - x1val
    dyVal = y2val - y1val
    dzVal = z2val - z1val

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    ld = sp.diff(lambHarmonic,lamb).subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

    return ld

def fixedEnergyCompute(ci,pos,posEq,r0val,Kval,lambdaVal):

    posp   = pos[ci]

    e = lambFixed.subs([(x1,posp[0]),(x2,posEq[0]),(y1,posp[1]),(y2,posEq[1]),(z1,posp[2]),(z2,posEq[2]),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

    return e

def fixedForceCompute(ci,pos,posEq,r0val,Kval,lambdaVal):

    posp   = pos[ci]

    fx = sp.diff(lambFixed,x1).subs([(x1,posp[0]),(x2,posEq[0]),(y1,posp[1]),(y2,posEq[1]),(z1,posp[2]),(z2,posEq[2]),(r0,r0val),(K,Kval),(lamb,lambdaVal)])
    fy = sp.diff(lambFixed,y1).subs([(x1,posp[0]),(x2,posEq[0]),(y1,posp[1]),(y2,posEq[1]),(z1,posp[2]),(z2,posEq[2]),(r0,r0val),(K,Kval),(lamb,lambdaVal)])
    fz = sp.diff(lambFixed,z1).subs([(x1,posp[0]),(x2,posEq[0]),(y1,posp[1]),(y2,posEq[1]),(z1,posp[2]),(z2,posEq[2]),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

    return np.array([float(fx),float(fy),float(fz)])

def fixedLambdaDerivativeCompute(ci,pos,posEq,r0val,Kval,lambdaVal):

    posp   = pos[ci]

    ld = sp.diff(lambFixed,lamb).subs([(x1,posp[0]),(x2,posEq[0]),(y1,posp[1]),(y2,posEq[1]),(z1,posp[2]),(z2,posEq[2]),(r0,r0val),(K,Kval),(lamb,lambdaVal)])

    return ld

def ljEnergyCompute(ids,pos,epsilonVal,sigmaVal,alphaVal,lambdaVal,cutOff=True):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    if(rval > 2.5*sigmaVal and cutOff):
        return 0.0

    e = lj.subs([(x1,x1val),(y1,y1val),(z1,z1val),
                 (x2,x2val),(y2,y2val),(z2,z2val),
                 (epsilon,epsilonVal),(sigma,sigmaVal),(alpha,alphaVal),(lamb,lambdaVal)])

    return e

def ljForceCompute(cid,ids,pos,epsilonVal,sigmaVal,alphaVal,lambdaVal,cutOff=True):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    dxVal = x2val - x1val
    dyVal = y2val - y1val
    dzVal = z2val - z1val

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    if(rval > 2.5*sigmaVal and cutOff):
        return np.array([0.0,0.0,0.0])

    fx = sp.diff(lj,x1).subs([(x1,x1val),(y1,y1val),(z1,z1val),
                              (x2,x2val),(y2,y2val),(z2,z2val),
                              (epsilon,epsilonVal),(sigma,sigmaVal),(alpha,alphaVal),(lamb,lambdaVal)])
    fy = sp.diff(lj,y1).subs([(x1,x1val),(y1,y1val),(z1,z1val),
                              (x2,x2val),(y2,y2val),(z2,z2val),
                              (epsilon,epsilonVal),(sigma,sigmaVal),(alpha,alphaVal),(lamb,lambdaVal)])
    fz = sp.diff(lj,z1).subs([(x1,x1val),(y1,y1val),(z1,z1val),
                              (x2,x2val),(y2,y2val),(z2,z2val),
                              (epsilon,epsilonVal),(sigma,sigmaVal),(alpha,alphaVal),(lamb,lambdaVal)])

    f = np.array([float(fx),float(fy),float(fz)])
    if cid == id2:
        return -f

    return f

def ljLambdaDerivativeCompute(cid,ids,pos,epsilonVal,sigmaVal,alphaVal,lambdaVal,cutOff=True):
    id1,id2 = ids
    pos1 = pos[id1]
    pos2 = pos[id2]

    x1val,y1val,z1val = pos1
    x2val,y2val,z2val = pos2

    dxVal = x2val - x1val
    dyVal = y2val - y1val
    dzVal = z2val - z1val

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    if(rval > 2.5*sigmaVal and cutOff):
        return 0.0

    rval = r.subs([(x1,x1val),(y1,y1val),(z1,z1val),(x2,x2val),(y2,y2val),(z2,z2val)])

    ld = sp.diff(lj,lamb).subs([(x1,x1val),(y1,y1val),(z1,z1val),
                                (x2,x2val),(y2,y2val),(z2,z2val),
                                (epsilon,epsilonVal),(sigma,sigmaVal),(alpha,alphaVal),(lamb,lambdaVal)])

    return ld

with open("./results/simulation.json") as f:
    sim = json.load(f)

ldTotalSim = 0.0
with open("./results/ti.dat") as f:
    for line in f:
        if "#" not in line:
            ldTotalSim += float(line)

###########################################

state = sim["state"]["data"]

ids = []
pos = []

for i,p in state:
    ids.append(i)
    pos.append(p)

ids = np.array(ids)
pos = np.array(pos)

###########################################

try:
    bondsData   = sim["topology"]["forceField"]["bonds"]["data"]
    ljBondsData = sim["topology"]["forceField"]["ljbonds"]["data"]
    fixedData   = sim["topology"]["forceField"]["fixed"]["data"]
except:
    bondsData   = []
    ljBondsData = []
    fixedData   = []

bonds   = []
ljBonds = []
fixed   = []

for b in bondsData:
    bonds.append([b[0], b[1]])

for b in ljBondsData:
    ljBonds.append([b[0], b[1]])

for f in fixedData:
    fixed.append([f[0], f[3]]) #id,poseq

try:
    r0val = bondsData[0][2]
    Kval  = bondsData[0][3]
    r0FixedVal = 0.0
    KFixedVal  = fixedData[0][2][0]
except:
    r0val = 0.0
    Kval  = 0.0
    r0FixedVal = 0.0
    KFixedVal  = 0.0

###########################################

potPath = "./results/pot.dat"
data = np.loadtxt(potPath,skiprows=2)

colBondsEnergy = 1
colBondsForces = 5

colFixedEnergy = 10
colFixedForces = 14

colLJEnergy = 19
colLJForces = 23

colLJBondsEnergy = 28
colLJBondsForces = 32

bondsEnergy = data[:,colBondsEnergy]
bondsForces = data[:,colBondsForces]

fixedEnergy = data[:,colFixedEnergy]
fixedForces = data[:,colFixedForces]

LJEnergy = data[:,colLJEnergy]
LJForces = data[:,colLJForces]

LJBondsEnergy = data[:,colLJBondsEnergy]
LJBondsForces = data[:,colLJBondsForces]

###########################################

#Compute energy bonds per particle

totalLambdaDerivative = 0.0

bondsEnergyPerParticle           = np.zeros(len(ids))
bondsForcesPerParticle           = np.zeros((len(ids),3))
bondsLambdaDerivativePerParticle = np.zeros(len(ids))
for i in ids:
    for b in bonds:
        if i in b:
            bondsEnergyPerParticle[i] += bondsEnergyCompute(b,pos,r0val,Kval,lambdaVal)/2.0

            f = bondsForceCompute(i,b,pos,r0val,Kval,lambdaVal)
            bondsForcesPerParticle[i] += f

            ld = bondsLambdaDerivativeCompute(i,b,pos,r0val,Kval,lambdaVal)/2.0
            bondsLambdaDerivativePerParticle[i] += ld
            totalLambdaDerivative += ld

bondsForcesMod = np.sqrt(bondsForcesPerParticle[:,0]**2 + bondsForcesPerParticle[:,1]**2 + bondsForcesPerParticle[:,2]**2)

ljBondsEnergyPerParticle = np.zeros(len(ids))
ljBondsForcesPerParticle = np.zeros((len(ids),3))
ljBondsLambdaDerivativePerParticle = np.zeros(len(ids))
for i in ids:
    for b in ljBonds:
        if i in b:
            j = b[0] if b[0] != i else b[1]
            ljBondsEnergyPerParticle[i] += ljEnergyCompute([i,j],pos,epsilonVal,r0val,alphaVal,lambdaVal,False)/2

            f = ljForceCompute(i,[i,j],pos,epsilonVal,r0val,alphaVal,lambdaVal,False)
            ljBondsForcesPerParticle[i] += f

            ld = ljLambdaDerivativeCompute(i,[i,j],pos,epsilonVal,r0val,alphaVal,lambdaVal,False)/2
            ljBondsLambdaDerivativePerParticle[i] += ld
            totalLambdaDerivative += ld

ljBondsForcesMod = np.sqrt(ljBondsForcesPerParticle[:,0]**2 + ljBondsForcesPerParticle[:,1]**2 + ljBondsForcesPerParticle[:,2]**2)

fixedEnergyPerParticle           = np.zeros(len(ids))
fixedForcesPerParticle           = np.zeros((len(ids),3))
fixedLambdaDerivativePerParticle = np.zeros(len(ids))
for i in ids:
    fixedEnergyPerParticle[i] = fixedEnergyCompute(i,pos,fixed[i][1],r0FixedVal,KFixedVal,lambdaVal)
    f = fixedForceCompute(i,pos,fixed[i][1],r0FixedVal,KFixedVal,lambdaVal)
    fixedForcesPerParticle[i] = f

    ld = fixedLambdaDerivativeCompute(i,pos,fixed[i][1],r0FixedVal,KFixedVal,lambdaVal)
    fixedLambdaDerivativePerParticle[i] = ld
    totalLambdaDerivative += ld

fixedForcesMod = np.sqrt(fixedForcesPerParticle[:,0]**2 + fixedForcesPerParticle[:,1]**2 + fixedForcesPerParticle[:,2]**2)

ljEnergyPerParticle           = np.zeros(len(ids))
ljForcesPerParticle           = np.zeros((len(ids),3))
ljLambdaDerivativePerParticle = np.zeros(len(ids))
for index,i in enumerate(ids):
    for j in ids[index+1:]:
        ljEnergyPerParticle[i] += ljEnergyCompute([i,j],pos,epsilonVal,sigmaVal,alphaVal,lambdaVal)/2
        ljEnergyPerParticle[j] += ljEnergyCompute([i,j],pos,epsilonVal,sigmaVal,alphaVal,lambdaVal)/2

        f = ljForceCompute(i,[i,j],pos,epsilonVal,sigmaVal,alphaVal,lambdaVal)
        ljForcesPerParticle[i] += f
        ljForcesPerParticle[j] -= f

        ld = ljLambdaDerivativeCompute(i,[i,j],pos,epsilonVal,sigmaVal,alphaVal,lambdaVal)
        ljLambdaDerivativePerParticle[i] += ld/2.0
        ljLambdaDerivativePerParticle[j] += ld/2.0
        totalLambdaDerivative += ld

ljForcesMod = np.sqrt(ljForcesPerParticle[:,0]**2 + ljForcesPerParticle[:,1]**2 + ljForcesPerParticle[:,2]**2)

#Compare bonds energy per particle
allCloseBE = True
for beTheo,beSim in zip(bondsEnergyPerParticle,bondsEnergy):
    if not np.isclose(beTheo,beSim):
        print("Error in bonds energy per particle (theo,sim):",beTheo,beSim)
        allCloseBE = False
        #break

allCloseBF = True
for bfTheo,bfSim in zip(bondsForcesMod,bondsForces):
    if not np.isclose(bfTheo,bfSim):
        print("Error in bonds forces per particle (theo,sim):",bfTheo,bfSim)
        allCloseBF = False
        #break

allCloseFE = True
for feTheo,feSim in zip(fixedEnergyPerParticle,fixedEnergy):
    if not np.isclose(feTheo,feSim):
        print("Error in fixed energy per particle (theo,sim):",feTheo,feSim)
        allCloseFE = False
        #break

allCloseFF = True
for ffTheo,ffSim in zip(fixedForcesMod,fixedForces):
    if not np.isclose(ffTheo,ffSim):
        print("Error in fixed forces per particle (theo,sim):",ffTheo,ffSim)
        allCloseFF = False
        #break

allCloseLJE = True
for ljeTheo,ljeSim in zip(ljEnergyPerParticle,LJEnergy):
    if not np.isclose(ljeTheo,ljeSim):
        print("Error in LJ energy per particle (theo,sim):",ljeTheo,ljeSim)
        allCloseLJE = False
        #break

allCloseLJF = True
for ljfTheo,ljfSim in zip(ljForcesMod,LJForces):
    if not np.isclose(ljfTheo,ljfSim):
        print("Error in LJ forces per particle (theo,sim):",ljfTheo,ljfSim)
        allCloseLJF = False
        #break

allCloseLJBE = True
for ljbeTheo,ljbeSim in zip(ljBondsEnergyPerParticle,LJBondsEnergy):
    if not np.isclose(ljbeTheo,ljbeSim):
        print("Error in LJ bonds energy per particle (theo,sim):",ljbeTheo,ljbeSim)
        allCloseLJBE = False
        #break

allCloseLJBF = True
for ljbfTheo,ljbfSim in zip(ljBondsForcesMod,LJBondsForces):
    if not np.isclose(ljbfTheo,ljbfSim):
        print("Error in LJ bonds forces per particle (theo,sim):",ljbfTheo,ljbfSim)
        allCloseLJBF = False
        #break

print("Does sim and theo match? ")
print("Bonds energy: ",allCloseBE)
print("Bonds forces: ",allCloseBF)
print("Fixed energy: ",allCloseFE)
print("Fixed forces: ",allCloseFF)
print("LJ energy: ",allCloseLJE)
print("LJ forces: ",allCloseLJF)
print("LJ bonds energy: ",allCloseLJBE)
print("LJ bonds forces: ",allCloseLJBF)

print("\nChecking lambda derivative...")
print("Total ld sim: ",ldTotalSim," Total ld theo: ",totalLambdaDerivative)
ldClose = np.isclose(float(ldTotalSim),float(totalLambdaDerivative))
print("Difference: ",ldTotalSim-totalLambdaDerivative,", close? ",ldClose)

allPassed = allCloseBE and allCloseBF and allCloseFE and allCloseFF and allCloseLJE and allCloseLJF and ldClose and allCloseLJBE and allCloseLJBF
print("\n\nTest passed? ",allPassed)
