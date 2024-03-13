import sys
import json

import numpy as np

from helixPotential import computePerParticleHelixEnergyForceTorque

from scipy.spatial.transform import Rotation

with open("pot.dat") as f:
    header = f.readline()
    _      = f.readline() # skip second line
    records = []
    for line in f:
        if line.startswith("#"):
            break
        records.append(line)

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

helixParams = sim["topology"]["forceField"]["helix"]["parameters"]

Eb = helixParams["E"]

Kb = helixParams["Kb"]
Ka = helixParams["Ka"]
Kd = helixParams["Kd"]

e_x = helixParams["e_x"]
e_y = helixParams["e_y"]
e_z = helixParams["e_z"]

R_H = np.asarray([e_x, e_y, e_z]).T

conn_next = helixParams["e_next"]
conn_prev = helixParams["e_prev"]

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
    e = float(info["helixEnergy"])
    p = np.asarray(info["position"])
    o = np.asarray(info["orientation"])
    q0,q1,q2,q3 = o
    o = Rotation.from_quat([q1,q2,q3,q0]).as_matrix()

    #print(o.T)
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

energyTheo,forcesTheo,torquesTheo = computePerParticleHelixEnergyForceTorque(positions,orientations,
                                                                             Eb,Kb,Ka,Kd,R_H,conn_next,conn_prev)

print("Total force: {}".format(totalForce))
print("Total torque: {}".format(totalTorque))

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
    print("All OK")
else:
    print("ERRORS, errorsCounter: {}".format(errorsCounter))

