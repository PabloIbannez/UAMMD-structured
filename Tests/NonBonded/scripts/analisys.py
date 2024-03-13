import sys,os

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

import json

#######################

COMPONENTS_PATH = os.path.join(os.getenv('UAMMD_PATH'),
                               'USCM/Components.json')

with open(COMPONENTS_PATH) as f:
    components = json.load(f)

#######################

if len(sys.argv) > 1:
    skip = int(sys.argv[1])
else:
    skip = 300

PARAMETERS_PATH = "./parameters.json"

with open(PARAMETERS_PATH, "r") as f:
    param = json.load(f)

nonBondedList = param["nonBondedList"]

units = param["units"]
if units == "none":
    kB      = 1.0
    ELECOEF = 1.0
elif units == "KcalMol_A":
    kB = 1.987191E-03
    ELECOEF = 332.0716

kBT      = kB*param["temperature"]
lambdaTI = param["lambda"]

print(f"Temperature: {param['temperature']} K")
print(f"Units: {units}")
print(f"Kb: {kB}")
print(f"ELECOEF: {ELECOEF}")
print(f"kBT: {kBT}")

chgProduct = param["charge"]*param["charge"]
radius     = param["radius"]

#######################

def Harmonic(r,K,r0):
    return 0.5*K*(r-r0)*(r-r0)

def MaxDistanceRestraint(r,K,maxDistance):
    return np.piecewise(r,[r<=maxDistance,r>maxDistance],[lambda r: 0.0,lambda r: Harmonic(r,K,maxDistance)])

def Atzberger(r,radius,epsilon,mu,chi,theta):
    return r*0.0

def DebyeHuckel(r,dielectricConstant,debyeLength,cutOffFactor):
    efactor = ELECOEF*chgProduct/dielectricConstant
    return efactor*np.exp(-r/debyeLength)/r

def DebyeHuckelSpheres(r,dielectricConstant,debyeLength,cutOffFactor):
    efactor = ELECOEF*chgProduct/dielectricConstant
    efactor = efactor/(1.0+radius/debyeLength)**2
    return efactor*np.exp(-(r-2.0*radius)/debyeLength)/r

def LennardJonesType1(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: 4.0*epsilon*(np.power(sigma/r,12)-np.power(sigma/r,6)),lambda r: 0.0])

def LennardJonesType2(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: epsilon*(np.power(sigma/r,12)-2.0*np.power(sigma/r,6)),lambda r: 0.0])

def LennardJonesType3(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: epsilon*(5.0*np.power(sigma/r,12)-6.0*np.power(sigma/r,10)),lambda r: 0.0])

def SplitLennardJones(r,epsilon_r,epsilon_a,epsilon,sigma,cutOffFactor):
    near_cut = np.power(2, 1/6)*sigma
    return np.where(r < near_cut,
                    LennardJonesType1(r,epsilon_r,sigma,cutOffFactor) + epsilon_r,
                    LennardJonesType1(r,epsilon_a*epsilon,sigma,cutOffFactor) + epsilon_a*epsilon)

def NonPolar(r,epsilon,sigma,cutOffFactor):
    if epsilon < 0:
        return LennardJonesType1(r,abs(epsilon),sigma,cutOffFactor)
    else:
        return np.piecewise(r,[r<=sigma*1.122462,r>sigma*1.122462],[lambda r: LennardJonesType1(r,epsilon,sigma,cutOffFactor)+ 2.0*epsilon,lambda r: -LennardJonesType1(r,epsilon,sigma,cutOffFactor)])

def DLVOType1(r,dielectricConstant,debyeLength,epsilon,sigma,cutOffNPFactor,cutOffDHFactor):
    return DebyeHuckelSpheres(r,dielectricConstant,debyeLength,cutOffDHFactor)+LennardJonesType1(r,epsilon,sigma,cutOffNPFactor)

def DLVOType2(r,dielectricConstant,debyeLength,epsilon,sigma,cutOffNPFactor,cutOffDHFactor):
    return DebyeHuckelSpheres(r,dielectricConstant,debyeLength,cutOffDHFactor)+LennardJonesType2(r,epsilon,sigma,cutOffNPFactor)

def DLVOType3(r,dielectricConstant,debyeLength,epsilon,sigma,cutOffNPFactor,cutOffDHFactor):
    return DebyeHuckelSpheres(r,dielectricConstant,debyeLength,cutOffDHFactor)+LennardJonesType3(r,epsilon,sigma,cutOffNPFactor)

def Clashed(r,gamma,lambda_):
    d = gamma*(2.0*radius)
    d2 = np.power(d,2)
    r2 = np.power(r,2)
    diff = d2-r2
    e = np.where(diff>0.0,diff,0.0)
    return lambda_*e

def KimHummer(r,dielectricConstant,debyeLength,epsilon,sigma,cutOffNPFactor,cutOffDHFactor,sasaModel):
    return NonPolar(r,epsilon,sigma,cutOffNPFactor)+DebyeHuckel(r,dielectricConstant,debyeLength,cutOffDHFactor)

def WCAType1(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma*1.122462,r>sigma*1.122462],[lambda r: LennardJonesType1(r,epsilon,sigma,cutOffFactor)+epsilon,lambda r: 0.0])

def WCAType2(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma,r>sigma],[lambda r: LennardJonesType2(r,epsilon,sigma,cutOffFactor)+epsilon,lambda r: 0.0])

def WCAType3(r,epsilon,sigma,cutOffFactor):
    return np.piecewise(r,[r<=sigma,r>sigma],[lambda r: LennardJonesType3(r,epsilon,sigma,cutOffFactor)+epsilon,lambda r: 0.0])

def GeneralLennardJonesType1(r,epsilon,sigma,cutOffFactor):
    if epsilon >= 0:
        return WCAType1(r,epsilon,sigma,cutOffFactor)
    else:
        return LennardJonesType1(r,np.abs(epsilon),sigma,cutOffFactor)

def GeneralLennardJonesType2(r,epsilon,sigma,cutOffFactor):
    if epsilon >= 0:
        return WCAType2(r,epsilon,sigma,cutOffFactor)
    else:
        return LennardJonesType2(r,np.abs(epsilon),sigma,cutOffFactor)

def GeneralLennardJonesType3(r,epsilon,sigma,cutOffFactor):
    if epsilon >= 0:
        return WCAType3(r,epsilon,sigma,cutOffFactor)
    else:
        return LennardJonesType3(r,np.abs(epsilon),sigma,cutOffFactor)

def Steric6(r,epsilon,sigma,cutOffFactor):
    return epsilon*(np.power(sigma/r,6))

def Steric12(r,epsilon,sigma,cutOffFactor):
    return epsilon*(np.power(sigma/r,12))

def Steric6SoftCore(r,epsilon,sigma,cutOffFactor,alpha):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: lambdaTI*epsilon/(alpha*(1.0-lambdaTI)**2+np.power(r/sigma,6)),lambda r: 0.0])

def Steric12SoftCore(r,epsilon,sigma,cutOffFactor,alpha):
    return np.piecewise(r,[r<=sigma*cutOffFactor,r>sigma*cutOffFactor],[lambda r: lambdaTI*epsilon/(alpha*(1.0-lambdaTI)**2+np.power(r/sigma,6))**2,lambda r: 0.0])

def LennardJonesSoftCoreType1(r,epsilon,sigma,cutOffFactor,alpha):
    return Steric12SoftCore(r,4.0*epsilon,sigma,cutOffFactor,alpha)-Steric6SoftCore(r,4.0*epsilon,sigma,cutOffFactor,alpha)

def LennardJonesSoftCoreType2(r,epsilon,sigma,cutOffFactor,alpha):
    return Steric12SoftCore(r,epsilon,sigma,cutOffFactor,alpha)-Steric6SoftCore(r,2.0*epsilon,sigma,cutOffFactor,alpha)

def boltz_bond2(x,func,funcParams):
    exponent =-func(x,**funcParams)/kBT
    return x*x*np.exp(exponent)

def boltz_bond2_bound(x,func,funcParams,boundK,boundMaxDistance):
    exponent =-func(x,**funcParams)/kBT
    exponent += -MaxDistanceRestraint(x,boundK,boundMaxDistance)/kBT
    return x*x*np.exp(exponent)

#######################

integrators = ["eulerMaruyama", "eulerMaruyamaRigid", "bbk", "gjf"]

for nonBondedCls,nonBondedSubCls,_ in components["Interactor"]["Pair"]:
    if nonBondedSubCls not in nonBondedList:
        print("[WARNING] NonBonded type {} not found".format(nonBondedSubCls))
        continue
    else:
        nonBondedInfo = nonBondedList[nonBondedSubCls]

    #Perform distance analysis
    for inteCount,inte in enumerate(integrators):

        print("\n\n#############################################\n\n")
        print(f"[INFO] Processing {nonBondedCls} {nonBondedSubCls} {inte}")

        file = f"./results/distances_{nonBondedSubCls}_{inte}.dat"
        potName = nonBondedSubCls

        print(f"Detected potential: {potName}")
        potFunc = eval(potName)

        #Bound

        boundParams = param["bond2bound"]

        boundType        = boundParams["type"]
        boundK           = boundParams["K"]
        boundMaxDistance = boundParams["distance"]

        if not os.path.isfile(file):
            print(f"[WARNING] File {file} not found")
            continue
        else:
            #Check if file is empty
            if os.stat(file).st_size == 0:
                print(f"[WARNING] File {file} is empty")
                continue
            print(f"[INFO] Analysing {file}")
            data = np.loadtxt(file, skiprows=skip)
            #If the size of data is zero, skip
            if data.size == 0:
                print(f"[WARNING] File {file} is empty")
                continue
            data = data.flatten()

            minDist = 0.0
            maxDist = boundMaxDistance+2.5
            r = np.linspace(minDist,maxDist,1000)

            typesParam = {}
            if "data" in nonBondedInfo:
                for index,lab in enumerate(nonBondedInfo["labels"]):
                    if "name" not in lab:
                        typesParam[lab] = nonBondedInfo["data"][0][index]

            if "lambda" in nonBondedInfo["parameters"].keys():
                nonBondedInfo["parameters"]["lambda_"] = nonBondedInfo["parameters"]["lambda"]
                del nonBondedInfo["parameters"]["lambda"]

            potParam = {**nonBondedInfo["parameters"], **typesParam}
            print(f"[INFO] Parameters: {potParam}")

            if(inteCount == 0):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(r,potFunc(r,**potParam)+MaxDistanceRestraint(r,boundK,boundMaxDistance),label="Potential")
                ax.set_ylim(-6.0,20.0)
                ax.set_xlabel("r")
                ax.set_ylabel("V(r)")
                #ax.set_title(f"Potential {nonBondedSubCls} {inte}, parameters: {potParam}")
                ax.set_title(f"Potential {nonBondedSubCls} {inte}")
                ax.legend()
                fig.savefig(f"./results/potential_{nonBondedSubCls}.png")

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(data, bins=200, density=True, label=f"{nonBondedSubCls} {inte}")

            if "sigma" in potParam.keys():
                limInf = potParam["sigma"]*0.5
            else:
                limInf = 0.0
            limSup = boundMaxDistance*2.0

            x = np.linspace(limInf,limSup,10000)
            prob = boltz_bond2_bound(x,potFunc,potParam,boundK,boundMaxDistance)

            #We have than prob=pot(x), we have to integrate it between limInf and limSup to compute Z
            Z = simps(prob,x)

            ax.plot(r,boltz_bond2_bound(r,potFunc,potParam,boundK,boundMaxDistance)/Z, label=f"{nonBondedSubCls} {inte}")

            ax.set_xlabel("r")
            ax.set_ylabel("P(r)")
            ax.set_title(f"Probability {nonBondedSubCls} {inte}")
            ax.legend()
            fig.savefig(f"./results/probability_{nonBondedSubCls}_{inte}.png")


