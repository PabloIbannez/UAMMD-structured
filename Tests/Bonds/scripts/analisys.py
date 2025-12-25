import sys,os

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simpson

from get_r0 import get_r0

import json

TORSION_IN_DEGREES = False

#######################

COMPONENTS_PATH = os.path.join(os.getenv('UAMMD_PATH'),'Components.json')

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

bondsList = param["bondsList"]

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

#######################

#Bond1

#Bond2

def Fene(r,K,r0,R0):
    R02 = R0*R0
    dr  = r-r0
    dr2 = dr*dr
    return -(K*R02/2.0)*np.log(1.0-dr2/R02)

def Harmonic(r,K,r0):
    return 0.5*K*(r-r0)*(r-r0)

def Morse(r,r0,E,D):
    dr = r-r0
    factor = 1.0-np.exp(-dr/D)
    return E*(factor*factor-1.0)

def MaxDistanceRestraint(r,K,maxDistance):
    return np.piecewise(r,[r<=maxDistance,r>maxDistance],[lambda r: 0.0,lambda r: Harmonic(r,K,maxDistance)])

def LennardJonesType1(r,epsilon,sigma):
    return 4.0*epsilon*(np.power(sigma/r,12)-np.power(sigma/r,6))

def LennardJonesType2(r,epsilon,sigma):
    return epsilon*(np.power(sigma/r,12)-2.0*np.power(sigma/r,6))

def LennardJonesType3(r,epsilon,sigma):
    return epsilon*(5.0*np.power(sigma/r,12)-6.0*np.power(sigma/r,10))

def WCAType2(r,epsilon,sigma):
    return np.piecewise(r,[r<=sigma,r>sigma],[lambda r: LennardJonesType2(r,epsilon,sigma)+epsilon,lambda r: 0.0])

def Gaussian(r,E,r0,D):
    dr = r-r0
    return -E*np.exp(-dr*dr/(2.0*D))

def LennardJonesGaussian(r,epsilon,sigma,D):
    return WCAType2(r,epsilon,sigma) + Gaussian(r,epsilon,sigma,D)

def LennardJonesKaranicolasBrooks(r,epsilon,sigma):
    return epsilon*(13.0*np.power(sigma/r,12)-18.0*np.power(sigma/r,10)+4.0*np.power(sigma/r,6))

def DebyeHuckel(r,chgProduct,dielectricConstant,debyeLength,cutOff):
    efactor = ELECOEF*chgProduct/dielectricConstant
    return np.piecewise(r,[r<=cutOff,r>cutOff],[lambda r: efactor*np.exp(-r/debyeLength)/r,lambda r: 0.0])

def MorseWCA(r,r0,E,D,eps0):
    return Morse(r,r0,E,D) + WCAType2(r,eps0,r0)

def Steric6(r,epsilon,sigma):
    return epsilon*(np.power(sigma/r,6))

def Steric12(r,epsilon,sigma):
    return epsilon*(np.power(sigma/r,12))

def LambdaHarmonic(r,K,r0):
    return 0.5*lambdaTI*K*(r-r0)*(r-r0)

def LennardJonesSoftCoreType1(r,epsilon,sigma,alpha):
    den = alpha*(1.0-lambdaTI)**2+np.power(r/sigma,6)
    e   = 4.0*lambdaTI*epsilon*(1.0/den**2-1.0/den)
    return e

def LennardJonesSoftCoreType2(r,epsilon,sigma,alpha):
    den = alpha*(1.0-lambdaTI)**2+np.power(r/sigma,6)
    e   = lambdaTI*epsilon*(1.0/den**2-2.0/den)
    return e

def boltz_bond2(x,func,funcParams):
    exponent =-func(x,**funcParams)/kBT
    return x*x*np.exp(exponent)

def boltz_bond2_bound(x,func,funcParams,boundK,boundMaxDistance):
    exponent =-func(x,**funcParams)/kBT
    exponent += -MaxDistanceRestraint(x,boundK,boundMaxDistance)/kBT
    return x*x*np.exp(exponent)

def BestChenHummerAngular(theta):
    #BCH params
    gamma = 0.1
    epsilon_alpha = 4.3

    theta_alpha = 1.6
    theta_beta  = 2.27

    k_alpha = 106.4
    k_beta  = 26.3

    adiff = theta-theta_alpha
    bdiff = theta-theta_beta
    adiff2 = adiff*adiff
    bdiff2 = bdiff*bdiff
    exp_alpha = np.exp(-gamma*(k_alpha*adiff2+epsilon_alpha))
    exp_beta  = np.exp(-gamma*(k_beta*bdiff2))
    expE = exp_alpha+exp_beta

    return -np.log(expE)/gamma

def KratkyPorod(theta,K):
    return K*(1.0+np.cos(theta))

def HarmonicAngular(theta,K,theta0):
    diff = theta-theta0
    return 0.5*K*diff*diff

def boltz_bond3(x,func,funcParams):
    exponent =-func(x,**funcParams)/kBT
    return np.sin(x)*np.exp(exponent)

def Dihedral(phi,n,K,phi0):
    diff = phi-phi0
    return K*(1.0+np.cos(n*diff))

def Dihedral(phi,K,n,phi0):
    return K*(1.0+np.cos(n*phi-phi0))

def Dihedral4(phi,K,phi0):
    K1,K2,K3,K4 = K
    phi0_1,phi0_2,phi0_3,phi0_4 = phi0
    e =K1*(1.0+np.cos(1*phi-phi0_1))
    e+=K2*(1.0+np.cos(2*phi-phi0_2))
    e+=K3*(1.0+np.cos(3*phi-phi0_3))
    e+=K4*(1.0+np.cos(4*phi-phi0_4))
    return e

def IDP_Fourier(phi):

    A = np.asarray([0.705,-0.313,-0.079,0.041])*0.593
    B = np.asarray([-0.175,-0.093,0.030,0.030])*0.593

    e =A[0]*np.cos(1*phi)+B[0]*np.sin(1*phi)
    e+=A[1]*np.cos(2*phi)+B[1]*np.sin(2*phi)
    e+=A[2]*np.cos(3*phi)+B[2]*np.sin(3*phi)
    e+=A[3]*np.cos(4*phi)+B[3]*np.sin(4*phi)

    return e

def boltz_bond4(x,func,funcParams):
    if TORSION_IN_DEGREES:
        exponent =-func(x*np.pi/180.0,**funcParams)/kBT
    else:
        exponent =-func(x,**funcParams)/kBT
    return np.exp(exponent)

#######################

integrators = ["eulerMaruyama", "eulerMaruyamaRigid", "bbk", "gjf"]

for bndCls,bndSubCls,_ in components["Interactor"]["Bonds"]:
    if bndSubCls not in bondsList:
        print("[WARNING] Bond {} not in bondsList".format(bndSubCls))
        continue
    else:
        bndInfo = bondsList[bndSubCls]

    if bndCls == "Bond1":
        continue
    elif bndCls == "Bond2":
        #Perform distance analysis
        for inteCount,inte in enumerate(integrators):

            print("\n\n#############################################\n\n")
            print(f"[INFO] Processing {bndCls} {bndSubCls} {inte}")

            file = f"./results/distances_{bndSubCls}_{inte}.dat"

            if "Common" in bndSubCls:
                #Get all text before "Common"
                potName = bndSubCls.split("Common")[0]
            else:
                potName = bndSubCls

            print(f"Detected potential: {potName}")
            potFunc = eval(potName)

            bound = False
            if bndInfo["needBounds"]:
                bound = True
                r0_bound = get_r0(bndCls,bndSubCls,bondsList)

                boundParams = param["bond2bound"]

                boundType   = boundParams["type"]
                boundK      = boundParams["K"]
                boundFactor = boundParams["factor"]

                boundMaxDistance = r0_bound*boundFactor

                print(f"Detected bound: {boundType}")
                print(f"r0 bound: {r0_bound}")
                print(f"Bound K: {boundK}")
                print(f"Bound maxDistance: {boundMaxDistance}")


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

                minDist = data.min()
                maxDist = data.max()
                r = np.linspace(minDist,maxDist,1000)

                #Merge dictionaries bndInfo["commonParam"] and bndInfo["bondParam"] int a single dictionary
                potParam = {**bndInfo["commonParams"], **bndInfo["bondParams"]}

                if(inteCount == 0):
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(r,potFunc(r,**potParam),label="Potential")
                    ax.set_xlabel("r")
                    ax.set_ylabel("V(r)")
                    ax.set_title(f"Potential {bndSubCls} {inte}, parameters: {potParam}")
                    ax.legend()
                    fig.savefig(f"./results/potential_{bndSubCls}.png")

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.hist(data, bins=200, density=True, label=f"{bndSubCls} {inte}")

                limInf = 0.0
                limSup = 100.0
                if "Fene" in potName:
                    r0 = potParam["r0"]
                    R0 = potParam["R0"]
                    limInf = r0-R0
                    limSup = r0+R0
                if "LennardJones" in potName or "WCA" in potName:
                    if "Morse" in potName:
                        sigma = potParam["r0"]
                    else:
                        sigma = potParam["sigma"]
                    limInf = sigma*0.5

                x = np.linspace(limInf,limSup,10000)

                if not bound:
                    prob = boltz_bond2(x,potFunc,potParam)
                else:
                    prob = boltz_bond2_bound(x,potFunc,potParam,boundK,boundMaxDistance)

                #We have than prob=pot(x), we have to integrate it between limInf and limSup to compute Z
                Z = simpson(prob,x)

                if not bound:
                    ax.plot(r,boltz_bond2(r,potFunc,potParam)/Z, label=f"{bndSubCls} {inte}")
                else:
                    ax.plot(r,boltz_bond2_bound(r,potFunc,potParam,boundK,boundMaxDistance)/Z, label=f"{bndSubCls} {inte}")

                ax.set_xlabel("r")
                ax.set_ylabel("P(r)")
                ax.set_title(f"Probability {bndSubCls} {inte}")
                ax.legend()
                fig.savefig(f"./results/probability_{bndSubCls}_{inte}.png")

    elif bndCls == "Bond3":
        #Perform distance analysis
        for inteCount,inte in enumerate(integrators):

            print("\n\n#############################################\n\n")
            print(f"[INFO] Processing {bndCls} {bndSubCls} {inte}")

            file = f"./results/angles_{bndSubCls}_{inte}.dat"

            if "Common" in bndSubCls:
                #Get all text before "Common"
                potName = bndSubCls.split("Common")[0]
            else:
                potName = bndSubCls

            print(f"Detected potential: {potName}")
            potFunc = eval(potName)

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
                if data.size == 0:
                    print(f"[WARNING] File {file} is empty")
                    continue
                data = data.flatten()

                minAngle = data.min()
                maxAngle = data.max()
                theta = np.linspace(minAngle,maxAngle,1000)

                #Merge dictionaries bndInfo["commonParam"] and bndInfo["bondParam"] int a single dictionary
                potParam = {**bndInfo["commonParams"], **bndInfo["bondParams"]}

                if(inteCount == 0):
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(theta,potFunc(theta,**potParam),label="Potential")
                    ax.set_xlabel("theta")
                    ax.set_ylabel("V(theta)")
                    ax.set_title(f"Potential {bndSubCls} {inte}, parameters: {potParam}")
                    ax.legend()
                    fig.savefig(f"./results/potential_{bndSubCls}.png")

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.hist(data, bins=200, density=True, label=f"{bndSubCls} {inte}")

                x = np.linspace(0,np.pi,10000)

                prob = boltz_bond3(x,potFunc,potParam)

                #We have than prob=pot(x), we have to integrate it between limInf and limSup to compute Z
                Z = simpson(prob,x)

                ax.plot(theta,boltz_bond3(theta,potFunc,potParam)/Z, label=f"{bndSubCls} {inte}")

                ax.set_xlabel("theta")
                ax.set_ylabel("P(theta)")
                ax.set_title(f"Probability {bndSubCls} {inte}")
                ax.legend()
                fig.savefig(f"./results/probability_{bndSubCls}_{inte}.png")

    elif bndCls == "Bond4":
        #Perform distance analysis
        for inteCount,inte in enumerate(integrators):

            print("\n\n#############################################\n\n")
            print(f"[INFO] Processing {bndCls} {bndSubCls} {inte}")

            file = f"./results/dihedrals_{bndSubCls}_{inte}.dat"

            if "Common" in bndSubCls:
                #Get all text before "Common"
                potName = bndSubCls.split("Common")[0]
            else:
                potName = bndSubCls

            print(f"Detected potential: {potName}")
            potFunc = eval(potName)

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
                if data.size == 0:
                    print(f"[WARNING] File {file} is empty")
                    continue
                if TORSION_IN_DEGREES:
                    data = np.rad2deg(data).flatten()
                else:
                    data = data.flatten()

                minPhi = data.min()
                maxPhi = data.max()
                phi = np.linspace(minPhi,maxPhi,1000)

                #Merge dictionaries bndInfo["commonParam"] and bndInfo["bondParam"] int a single dictionary
                potParam = {**bndInfo["commonParams"], **bndInfo["bondParams"]}

                if(inteCount == 0):
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(phi,potFunc(phi,**potParam),label="Potential")
                    ax.set_xlabel("phi")
                    ax.set_ylabel("V(phi)")
                    ax.set_title(f"Potential {bndSubCls} {inte}, parameters: {potParam}")
                    ax.legend()
                    fig.savefig(f"./results/potential_{bndSubCls}.png")

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.hist(data, bins=200, density=True, label=f"{bndSubCls} {inte}")

                if TORSION_IN_DEGREES:
                    x = np.linspace(-180,180,10000)
                else:
                    x = np.linspace(-np.pi,np.pi,10000)

                prob = boltz_bond4(x,potFunc,potParam)

                #We have than prob=pot(x), we have to integrate it between limInf and limSup to compute Z
                Z = simpson(prob,x)

                ax.plot(phi,boltz_bond4(phi,potFunc,potParam)/Z, label=f"{bndSubCls} {inte}")

                ax.set_xlabel("phi")
                ax.set_ylabel("P(phi)")
                ax.set_title(f"Probability {bndSubCls} {inte}")
                ax.legend()
                fig.savefig(f"./results/probability_{bndSubCls}_{inte}.png")

    else:
        print(f"[ERROR] Bond class {bndCls} not supported")
        sys.exit(1)

