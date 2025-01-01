import os,sys
import json

import numpy as np
import quaternion as qt

from potential import computePerParticleEnergyForceTorque

sys.path.append(os.path.join(sys.path[0],'..','..','scripts'))
from utils import getVx,getVy,getVz

with open("pot.dat") as f:
    header = f.readline()
    _      = f.readline() # skip second line
    records  = []
    records1 = []
    frame = 0
    for line in f:
        if line.startswith("#"):
            frame += 1
            continue

        if frame == 0:
            records.append(line)
        elif frame == 1:
            records1.append(line)
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

partInfo1 = {}
for r in records1:
    r = r.split()
    partInfo1[int(r[0])] = {e:float(v) for e,v in zip(entries[1:], r[1:])}

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

cosrapParamsLabels = sim["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["labels"]
cosrapParamsData   = sim["topology"]["forceField"]["cosrap"]["patchesTopology"]["forceField"]["cosrap"]["data"][0]

cosrapParams = {l:d for l,d in zip(cosrapParamsLabels, cosrapParamsData)}

Eb = cosrapParams["E"]

Kswt = cosrapParams["Kswt"]
Krap = cosrapParams["Krap"]
rc   = cosrapParams["rc"]

Rq = cosrapParams['R']
R  = qt.quaternion(*Rq)
R  = qt.as_rotation_matrix(R)

print("R input:\n{}".format(R))

rxx = R[0,0]
rxy = R[1,0]
rxz = R[2,0]

ryx = R[0,1]
ryy = R[1,1]
ryz = R[2,1]

rzx = R[0,2]
rzy = R[1,2]
rzz = R[2,2]

connections = sim["topology"]["forceField"]["cosrap"]["patchesState"]["data"]

conn_S = connections[0][1]
conn_E = connections[1][1]

print("Connection start: {}".format(conn_S))
print("Connection end:   {}".format(conn_E))

state = sim["state"]["data"]

for s in state:
    id_,pos,ori = s
    partInfo[id_]["position"]    = pos
    partInfo[id_]["orientation"] = ori

    partInfo1[id_]["position"]    = pos
    partInfo1[id_]["orientation"] = ori


# Frame 0

totalEnergy = 0.0
totalForce  = np.asarray([0.,0.,0.])
totalTorque = np.asarray([0.,0.,0.])

energy       = []
positions    = []
orientations = []
forces       = []
torques      = []

for id_,info in partInfo.items():
    e = float(info["cosrapEnergy"])
    p = np.asarray(info["position"])
    o = np.asarray(info["orientation"])

    vx = getVx(o)
    vy = getVy(o)
    vz = getVz(o)

    o = np.asarray([vx,vy,vz]).T

    f = np.asarray([info["cosrapForceX"], info["cosrapForceY"], info["cosrapForceZ"]])
    t = np.asarray([info["cosrapTorqueX"], info["cosrapTorqueY"], info["cosrapTorqueZ"]])
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

energyTheo,forcesTheo,torquesTheo,bonds = computePerParticleEnergyForceTorque(positions,orientations,
                                                                              None,
                                                                              box,
                                                                              Eb,
                                                                              sigma,
                                                                              Kswt,Krap,
                                                                              rc,
                                                                              conn_S,conn_E,
                                                                              rxx,rxy,rxz,
                                                                              ryx,ryy,ryz,
                                                                              rzx,rzy,rzz)

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
    print("Frame0: All OK")
else:
    print("Frame0: ERRORS, errorsCounter: {}".format(errorsCounter))

# Frame 1
print("################################################")
print("################################################")
print("################################################")
print("################################################")

totalEnergy = 0.0
totalForce  = np.asarray([0.,0.,0.])
totalTorque = np.asarray([0.,0.,0.])

energy       = []
positions    = []
orientations = []
forces       = []
torques      = []

for id_,info in partInfo1.items():
    e = float(info["cosrapEnergy"])
    p = np.asarray(info["position"])

    o = np.asarray(info["orientation"])

    vx = getVx(o)
    vy = getVy(o)
    vz = getVz(o)

    o = np.asarray([vx,vy,vz]).T

    f = np.asarray([info["cosrapForceX"], info["cosrapForceY"], info["cosrapForceZ"]])
    t = np.asarray([info["cosrapTorqueX"], info["cosrapTorqueY"], info["cosrapTorqueZ"]])
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

print("Total force:  {}".format(totalForce))
print("Total torque: {}".format(totalTorque))

energyTheo,forcesTheo,torquesTheo,_ = computePerParticleEnergyForceTorque(positions,orientations,
                                                                          bonds,
                                                                          box,
                                                                          Eb,
                                                                          sigma,
                                                                          Kswt,Krap,
                                                                          rc,
                                                                          conn_S,conn_E,
                                                                          rxx,rxy,rxz,
                                                                          ryx,ryy,ryz,
                                                                          rzx,rzy,rzz)

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
    print("Frame1: All OK")
else:
    print("Frame1: ERRORS, errorsCounter: {}".format(errorsCounter))

    #for i,b in enumerate(bonds):
    #    if b < i or b == -1:
    #        continue
    #    print(i,b)

