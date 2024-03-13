import os,sys

import numpy as np
import matplotlib.pyplot as plt

import sympy as sp
# Import lambdify
from sympy.utilities.lambdify import lambdify

# Import kdtree from scipy
from MDAnalysis.lib.pkdtree import PeriodicKDTree

from collections import OrderedDict

from tqdm import tqdm

Eb_sym = sp.symbols("Eb_sym")

r_s_sym = sp.symbols("r_s_sym")
rc_sym = sp.symbols("rc_sym")

# Define position vectors for the start and end particles
r_start = sp.Matrix(sp.symbols('r_sx r_sy r_sz'))
r_end   = sp.Matrix(sp.symbols('r_ex r_ey r_ez'))

# Define orientation vectors for the start particle
s_x = sp.Matrix(sp.symbols('s_xx s_xy s_xz'))
s_y = sp.Matrix(sp.symbols('s_yx s_yy s_yz'))
s_z = sp.Matrix(sp.symbols('s_zx s_zy s_zz'))

# Define orientation vectors for the end particle
e_x = sp.Matrix(sp.symbols('e_xx e_xy e_xz'))
e_y = sp.Matrix(sp.symbols('e_yx e_yy e_yz'))
e_z = sp.Matrix(sp.symbols('e_zx e_zy e_zz'))

#
c_i = sp.Matrix(sp.symbols('c_ix c_iy c_iz'))
c_j = sp.Matrix(sp.symbols('c_jx c_jy c_jz'))

def analyticPotential():

    #####################################################################

    # Lab frame connections c_i_lab = S@c_i
    c_i_lab_x = s_x[0]*c_i[0] + s_y[0]*c_i[1] + s_z[0]*c_i[2]
    c_i_lab_y = s_x[1]*c_i[0] + s_y[1]*c_i[1] + s_z[1]*c_i[2]
    c_i_lab_z = s_x[2]*c_i[0] + s_y[2]*c_i[1] + s_z[2]*c_i[2]

    c_i_lab = sp.Matrix([c_i_lab_x,c_i_lab_y,c_i_lab_z])

    # Lab frame connections c_j_lab = E@c_j
    c_j_lab_x = e_x[0]*c_j[0] + e_y[0]*c_j[1] + e_z[0]*c_j[2]
    c_j_lab_y = e_x[1]*c_j[0] + e_y[1]*c_j[1] + e_z[1]*c_j[2]
    c_j_lab_z = e_x[2]*c_j[0] + e_y[2]*c_j[1] + e_z[2]*c_j[2]

    c_j_lab = sp.Matrix([c_j_lab_x,c_j_lab_y,c_j_lab_z])

    c_i_lab = c_i_lab + r_start
    c_j_lab = c_j_lab + r_end

    # Compute distance between connections

    r = sp.sqrt((c_i_lab[0]-c_j_lab[0])**2 + (c_i_lab[1]-c_j_lab[1])**2 + (c_i_lab[2]-c_j_lab[2])**2)

    normalized_r = (r - r_s_sym)/(rc_sym - r_s_sym)

    Ub = sp.Piecewise((-1.0,r < r_s_sym),
                      (-sp.cos(sp.pi*normalized_r),(r >= r_s_sym) & (r <= rc_sym)),
                      ( 1.0,r > rc_sym))

    Ub = (1.0 - Ub)/2.0

    ####################################################

    U = -Eb_sym*Ub

    #print("Computing S force ...")
    Uforce   = -U.diff(r_start)
    #print("Computing S torque ...")
    Utorque  = - s_x.cross(U.diff(s_x)) - s_y.cross(U.diff(s_y)) - s_z.cross(U.diff(s_z))

    ##################################

    U_dict = OrderedDict()

    U_dict['r'] = r
    U_dict['U'] = U
    U_dict['Uforce'] = Uforce
    U_dict['Utorque'] = Utorque

    # Lambdify all the expressions
    var_list = [Eb_sym,
                r_s_sym,
                rc_sym,
                r_start[0],r_start[1],r_start[2],
                r_end[0],r_end[1],r_end[2],
                s_x[0],s_x[1],s_x[2],
                s_y[0],s_y[1],s_y[2],
                s_z[0],s_z[1],s_z[2],
                e_x[0],e_x[1],e_x[2],
                e_y[0],e_y[1],e_y[2],
                e_z[0],e_z[1],e_z[2],
                c_i[0],c_i[1],c_i[2],
                c_j[0],c_j[1],c_j[2]]

    print("Lambdifying expressions...")
    for key in U_dict:
        U_dict[key] = lambdify(var_list,U_dict[key],modules='numpy')
    print("Done.")

    return U_dict


def evalAnalytic(analyticExpression,
                 ti,tj,
                 cn_i,cn_j,
                 pos_i,pos_j, # Parent
                 ori_i,ori_j,Eb,
                 r_s,rc):

    cn_s  = cn_i
    pos_s = pos_i
    ori_s = ori_i

    cn_e  = cn_j
    pos_e = pos_j
    ori_e = ori_j

    analyticEval = analyticExpression(Eb,
                                      r_s,rc,
                                      pos_s[0],pos_s[1],pos_s[2],
                                      pos_e[0],pos_e[1],pos_e[2],
                                      ori_s[0,0],ori_s[1,0],ori_s[2,0],
                                      ori_s[0,1],ori_s[1,1],ori_s[2,1],
                                      ori_s[0,2],ori_s[1,2],ori_s[2,2],
                                      ori_e[0,0],ori_e[1,0],ori_e[2,0],
                                      ori_e[0,1],ori_e[1,1],ori_e[2,1],
                                      ori_e[0,2],ori_e[1,2],ori_e[2,2],
                                      cn_s[0],cn_s[1],cn_s[2],
                                      cn_e[0],cn_e[1],cn_e[2])

    return analyticEval

def computePerParticleEnergyForceTorque(pos,ori,
                                        bondList,
                                        box,
                                        Eb,
                                        sigma,
                                        r_s,rc,
                                        conn_X_up,conn_X_down,
                                        conn_Y_up,conn_Y_down,
                                        conn_Z_up,conn_Z_down):


    U_dict = analyticPotential()
    N = pos.shape[0]

    energy  = []
    forces  = []
    torques = []

    for i in range(N):
        energy.append(0.0)
        forces.append(np.zeros(3))
        torques.append(np.zeros(3))

    ## Create patches

    patches    = []
    patchesTpy = []
    parentId   = []

    tpy2conn = {}

    tpy2conn["XU"] = conn_X_up
    tpy2conn["XD"] = conn_X_down
    tpy2conn["YU"] = conn_Y_up
    tpy2conn["YD"] = conn_Y_down
    tpy2conn["ZU"] = conn_Z_up
    tpy2conn["ZD"] = conn_Z_down

    for i in range(N):
        p = pos[i]
        o = ori[i]

        c_Xu = o@conn_X_up + p
        c_Xd = o@conn_X_down + p
        c_Yu = o@conn_Y_up + p
        c_Yd = o@conn_Y_down + p
        c_Zu = o@conn_Z_up + p
        c_Zd = o@conn_Z_down + p

        patches.append(c_Xu)
        patches.append(c_Xd)
        patches.append(c_Yu)
        patches.append(c_Yd)
        patches.append(c_Zu)
        patches.append(c_Zd)

        patchesTpy.append("XU")
        patchesTpy.append("XD")
        patchesTpy.append("YU")
        patchesTpy.append("YD")
        patchesTpy.append("ZU")
        patchesTpy.append("ZD")

        parentId.append(i)
        parentId.append(i)
        parentId.append(i)
        parentId.append(i)
        parentId.append(i)
        parentId.append(i)

    patch2parent = {}
    for i in range(len(patches)):
        patch2parent[i] = parentId[i]

    # Create KDtree using patches and rc
    tree = PeriodicKDTree(box=np.asarray(box+[90,90,90],dtype=np.float32))
    tree.set_coords(np.asarray(patches,dtype=np.float32),cutoff=rc)

    # Query pairs
    pairs = tree.search_pairs(rc)

    # Convert pairs to a list of tuples
    pairs = [(p[0],p[1]) for p in pairs]

    if bondList is not None:
        pairsTmp = []
        isBonded = []
        # Check all bonds in bondList are in pairs
        for i,b in enumerate(bondList):
            #print(i,b)
            if b < i or b == -1:
                continue
            p = (i,b)
            pairsTmp.append(p)
            isBonded.append(i)
            isBonded.append(b)

            if p not in pairs:
                raise RuntimeError("Bond "+str(p)+" not in pairs")

        if len(isBonded) != len(set(isBonded)):
            raise RuntimeError("Bond list contains duplicate bonds")

        for i,j in pairs:
            if i not in isBonded and j not in isBonded:
                #print("Adding pair: ",i,j)
                pairsTmp.append((i,j))

        for i,j in pairsTmp:
            if i > j:
                raise RuntimeError("i > j")

        pairs = pairsTmp.copy()

    patch2neigh = {}
    for p in range(len(patches)):
        patch2neigh[p] = []

    for p in pairs:
        patch2neigh[p[0]].append(p[1])
        patch2neigh[p[1]].append(p[0])

    bonds = []
    for pi in patch2neigh.keys():
        bonds.append(-1)
    bondsEnergy = []
    for pi in patch2neigh.keys():
        bondsEnergy.append(0.0)

    #logFile = open("dataTheo.dat","w")
    for pi in patch2neigh.keys():

        parentId = patch2parent[pi]

        neigh = patch2neigh[pi]

        pos_i = pos[patch2parent[pi]]
        ori_i = ori[patch2parent[pi]]

        tpy_i = patchesTpy[pi]

        cn_i = tpy2conn[tpy_i]

        localEnergy = 0.0

        localForce  = np.zeros(3)
        localTorque = np.zeros(3)

        minimalEnergy = 0.0

        for pj in neigh:

            pos_j = pos[patch2parent[pj]]
            ori_j = ori[patch2parent[pj]]

            tpy_j = patchesTpy[pj]

            cn_j = tpy2conn[tpy_j]

            bondEnergy = evalAnalytic(U_dict['U'],
                                      tpy_i,tpy_j,
                                      cn_i,cn_j,
                                      pos_i,pos_j,
                                      ori_i,ori_j,Eb,
                                      r_s,rc)

            if bondEnergy/2.0 < minimalEnergy:
                minimalEnergy = bondEnergy/2.0
                bonds[pi] = pj
                bondsEnergy[pi] = bondEnergy/2.0


            localEnergy+= bondEnergy/2

            fs = evalAnalytic(U_dict['Uforce'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,Eb,
                              r_s,rc)

            localForce  += np.asarray([float(x) for x in fs])

            ts = evalAnalytic(U_dict['Utorque'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,Eb,
                              r_s,rc)

            localTorque += np.asarray([float(x) for x in ts])

        energy[parentId]  += localEnergy
        forces[parentId]  += localForce
        torques[parentId] += localTorque

    #logFile.close()

    for pi in patch2neigh.keys():
        if bondsEnergy[pi] < 0.0:
            if (pi == bonds[bonds[pi]]):
                pass
            else:
                bonds[pi] = -1
        else:
            bonds[pi] = -1

    # Check that if pi is bonded with pj then pj is bonded with pi
    allOk = True
    for pi in patch2neigh.keys():
        pj = bonds[pi]
        if pj != -1:
            if pi != bonds[pj]:
                allOk = False
                print(f"Error: {pi} ({patch2parent[pi]}) is bonded with {pj} ({patch2parent[pj]}) but {pj} ({patch2parent[pj]}) is bonded with {bonds[pj]} ({patch2parent[bonds[pj]]})")

    if(allOk == False):
        raise ValueError("Error: Bonded pair not symmetric")

    return energy,forces,torques,bonds





