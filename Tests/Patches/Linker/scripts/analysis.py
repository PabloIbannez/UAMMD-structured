import os,sys
import json

import numpy as np

from linkerPotential import computePerParticleLinkerEnergyForceTorque

def getVx(quat):
    q0 = quat[0]
    q1 = quat[1]
    q2 = quat[2]
    q3 = quat[3]

    return [q0*q0+q1*q1-q2*q2-q3*q3,2.0*(q1*q2+q0*q3),2.0*(q1*q3-q0*q2)]

def getVy(quat):
    q0 = quat[0]
    q1 = quat[1]
    q2 = quat[2]
    q3 = quat[3]

    return [2.0*(q1*q2-q0*q3),q0*q0-q1*q1+q2*q2-q3*q3,2.0*(q2*q3+q0*q1)]

def getVz(quat):
    q0 = quat[0]
    q1 = quat[1]
    q2 = quat[2]
    q3 = quat[3]

    return [2.0*(q1*q3+q0*q2),2.0*(q2*q3-q0*q1),q0*q0-q1*q1-q2*q2+q3*q3]

with open("pot.dat") as f:
    header = f.readline()
    _      = f.readline() # skip second line
    records  = []
    frame = 0
    for line in f:
        if line.startswith("#"):
            frame += 1
            continue

        if frame == 0:
            records.append(line)
        else:
            break

entries = [h for h in header.split()[1:]]

#Check if the first entry is the id
if entries[0] != "id":
    print("The first entry of the header should be 'id'. But it is: {}".format(entries[0]))
    sys.exit(1)

partInfo = {}
for r in records:
    r = r.split()
    partInfo[int(r[0])] = {e:float(v) for e,v in zip(entries[1:], r[1:])}

with open("simulation.json") as f:
    sim = json.load(f)

ensembleLabels = sim["global"]["ensemble"]["labels"]
ensembleData   = sim["global"]["ensemble"]["data"][0]

ensemble = {l:v for l,v in zip(ensembleLabels, ensembleData)}

box = ensemble["box"]

print("Box: {}".format(box))

monomerTypeLabels = sim["global"]["types"]["labels"]
monomerTypeData   = sim["global"]["types"]["data"][0]

monomerType = {l:d for l,d in zip(monomerTypeLabels, monomerTypeData)}

sigma = monomerType["radius"]*2.0

linkerParamsLabels = sim["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["labels"]
linkerParamsData   = sim["topology"]["forceField"]["linker"]["patchesTopology"]["forceField"]["linker"]["data"][0]

linkerParams = {l:d for l,d in zip(linkerParamsLabels, linkerParamsData)}

E  = linkerParams["epsilon"]
D  = linkerParams["sigma"]

linkerPos = sim["topology"]["forceField"]["linker"]["patchesState"]["data"][0][1]

state = sim["state"]["data"]

for s in state:
    id_,pos,ori = s
    partInfo[id_]["position"]    = pos
    partInfo[id_]["orientation"] = ori

totalEnergy = 0.0
totalForce  = np.asarray([0.,0.,0.])
totalTorque = np.asarray([0.,0.,0.])

energy       = []
positions    = []
orientations = []
forces       = []
torques      = []

for id_,info in partInfo.items():
    e = float(info["linkerEnergy"])
    p = np.asarray(info["position"])
    o = np.asarray(info["orientation"])

    vx = getVx(o)
    vy = getVy(o)
    vz = getVz(o)

    o = np.asarray([vx,vy,vz]).T

    f = np.asarray([info["linkerForceX"], info["linkerForceY"], info["linkerForceZ"]])
    t = np.asarray([info["linkerTorqueX"], info["linkerTorqueY"], info["linkerTorqueZ"]])
    totalEnergy += e
    totalForce  += f
    totalTorque += t + np.cross(p,f)

    energy.append(e)
    positions.append(p)
    orientations.append(o)
    forces.append(f)
    torques.append(t)

positions    = np.asarray(positions)
orientations = np.asarray(orientations)

print("Total force: {}".format(totalForce))
print("Total torque: {}".format(totalTorque))

energyTheo,forcesTheo,torquesTheo = computePerParticleLinkerEnergyForceTorque(positions,orientations,
                                                                                    E,D,linkerPos,0.0)

N = len(positions)

allOk = True
errorsCounter = 0

print("Energy")
for i in range(N):
    e     = energy[i]
    eTheo = energyTheo[i]
    diff  = np.abs(e - eTheo)
    equal = diff < 2e-3
    equal = "OK" if equal else "ERROR"
    if equal != "OK":
        allOk = False
        errorsCounter += 1
        print(f"Id {i}, energy: {e}, energyTheo: {eTheo}, equal: {equal}, maxDiff: {diff}")
    else:
        print(f"Id {i}, energy: {e}, energyTheo: {eTheo}, equal: {equal}")

print("Force")
for i in range(N):
    f     = forces[i]
    fTheo = forcesTheo[i]
    diff  = np.abs(f - fTheo)
    equal = (diff < 2e-3).all()
    equal = "OK" if equal else "ERROR"
    if equal != "OK":
        allOk = False
        errorsCounter += 1
        print(f"Id {i}, force: {f}, forceTheo: {fTheo}, equal: {equal}, maxDiff: {np.max(diff)}")
    else:
        print(f"Id {i}, force: {f}, forceTheo: {fTheo}, equal: {equal}")

print("Torque")
for i in range(N):
    t     = torques[i]
    tTheo = torquesTheo[i]
    diff  = np.abs(t - tTheo)
    equal = (diff < 2e-3).all()
    equal = "OK" if equal else "ERROR"
    if equal != "OK":
        allOk = False
        errorsCounter += 1
        print(f"Id {i}, torque: {t}, torqueTheo: {tTheo}, equal: {equal}, maxDiff: {np.max(diff)}")
    else:
        print(f"Id {i}, torque: {t}, torqueTheo: {tTheo}, equal: {equal}")

if allOk:
    print("Frame: All OK")
else:
    print("Frame: ERRORS, errorsCounter: {}".format(errorsCounter))

