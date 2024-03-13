import os,sys
import json

import numpy as np

from helixPotential import computePerParticleHelixEnergyForceTorque

sys.path.append(os.path.join(sys.path[0],'..','..','scripts'))
from helixGeneration import getVx,getVy,getVz

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

helixParamsLabels = sim["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["labels"]
helixParamsData   = sim["topology"]["forceField"]["helix"]["patchesTopology"]["forceField"]["helix"]["data"][0]

helixParams = {l:d for l,d in zip(helixParamsLabels, helixParamsData)}

Eb = helixParams["Eb"]

Kb = helixParams["Kb"]
Ka = helixParams["Ka"]
Kd = helixParams["Kd"]

rc = helixParams["rc"]

e_x = helixParams["e_x"]
e_y = helixParams["e_y"]
e_z = helixParams["e_z"]

R_H = np.asarray([e_x, e_y, e_z]).T

connections = sim["topology"]["forceField"]["helix"]["patchesState"]["data"]

conn_prev = connections[0][1]
conn_next = connections[1][1]

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
    e = float(info["helixEnergy"])
    p = np.asarray(info["position"])
    o = np.asarray(info["orientation"])

    vx = getVx(o)
    vy = getVy(o)
    vz = getVz(o)

    o = np.asarray([vx,vy,vz]).T

    f = np.asarray([info["helixForceX"], info["helixForceY"], info["helixForceZ"]])
    t = np.asarray([info["helixTorqueX"], info["helixTorqueY"], info["helixTorqueZ"]])
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

energyTheo,forcesTheo,torquesTheo,bonds = computePerParticleHelixEnergyForceTorque(positions,orientations,
                                                                                   None,
                                                                                   box,
                                                                                   Eb,
                                                                                   sigma,rc,
                                                                                   Kb,Ka,Kd,R_H,conn_next,conn_prev)

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
    e = float(info["helixEnergy"])
    p = np.asarray(info["position"])

    o = np.asarray(info["orientation"])

    vx = getVx(o)
    vy = getVy(o)
    vz = getVz(o)

    o = np.asarray([vx,vy,vz]).T

    f = np.asarray([info["helixForceX"], info["helixForceY"], info["helixForceZ"]])
    t = np.asarray([info["helixTorqueX"], info["helixTorqueY"], info["helixTorqueZ"]])
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

energyTheo,forcesTheo,torquesTheo,_ = computePerParticleHelixEnergyForceTorque(positions,orientations,
                                                                               bonds,
                                                                               box,
                                                                               Eb,
                                                                               sigma,rc,
                                                                               Kb,Ka,Kd,R_H,conn_next,conn_prev)

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

