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

Kswt_sym  = sp.symbols("Kswt_sym")
Krap_sym  = sp.symbols("Krap_sym")
rc_sym    = sp.symbols("rc_sym")

# Define the rotation matrix
Rxx, Rxy, Rxz = sp.symbols('R_xx R_xy R_xz')
Ryx, Ryy, Ryz = sp.symbols('R_yx R_yy R_yz')
Rzx, Rzy, Rzz = sp.symbols('R_zx R_zy R_zz')

R = sp.Matrix([[Rxx, Ryx, Rzx],
               [Rxy, Ryy, Rzy],
               [Rxz, Ryz, Rzz]])

# Define position vectors for the start and end particles
r_start = sp.Matrix(sp.symbols('r_sx r_sy r_sz'))
r_end   = sp.Matrix(sp.symbols('r_ex r_ey r_ez'))

# Define orientation vectors for the start particle
Sxx, Sxy, Sxz = sp.symbols('S_xx S_xy S_xz')
Syx, Syy, Syz = sp.symbols('S_yx S_yy S_yz')
Szx, Szy, Szz = sp.symbols('S_zx S_zy S_zz')

s_x = sp.Matrix([Sxx, Sxy, Sxz])
s_y = sp.Matrix([Syx, Syy, Syz])
s_z = sp.Matrix([Szx, Szy, Szz])

S = sp.Matrix([[Sxx, Syx, Szx],
               [Sxy, Syy, Szy],
               [Sxz, Syz, Szz]])

# Define orientation vectors for the end particle
Exx, Exy, Exz = sp.symbols('E_xx E_xy E_xz')
Eyx, Eyy, Eyz = sp.symbols('E_yx E_yy E_yz')
Ezx, Ezy, Ezz = sp.symbols('E_zx E_zy E_zz')

E = sp.Matrix([[Exx, Eyx, Ezx],
               [Exy, Eyy, Ezy],
               [Exz, Eyz, Ezz]])

#
c_i = sp.Matrix(sp.symbols('c_ix c_iy c_iz'))
c_j = sp.Matrix(sp.symbols('c_jx c_jy c_jz'))

def analyticPotential():

    #####################################################################

    # Lab frame connections c_i_lab = S@c_i
    c_i_lab = S@c_i

    # Lab frame connections c_j_lab = E@c_j
    c_j_lab = E@c_j

    c_i_lab = c_i_lab + r_start
    c_j_lab = c_j_lab + r_end

    # Compute distance between connections

    r = sp.sqrt((c_i_lab[0]-c_j_lab[0])**2 + (c_i_lab[1]-c_j_lab[1])**2 + (c_i_lab[2]-c_j_lab[2])**2)

    normalized_r = r/rc_sym

    swt = (1.0/2.0)*(1.0 + sp.cos(sp.pi*normalized_r))

    def g(x,K):
        return x/(x + K*(1.0 - x))

    Ub = sp.Piecewise((g(swt,Kswt_sym),r < rc_sym),
                      ( 0.0,r > rc_sym))

    ####################################################

    frob = (E.T@S@R).trace()
    frob = (1.0 + frob)/4.0

    Urap  = g(frob,Krap_sym)

    ####################################################

    U = -Eb_sym*Ub*Urap

    #print("Computing S force ...")
    Uforce   = -U.diff(r_start)
    #print("Computing S torque ...")
    Utorque  = - s_x.cross(U.diff(s_x)) - s_y.cross(U.diff(s_y)) - s_z.cross(U.diff(s_z))

    ##################################

    U_dict = OrderedDict()

    U_dict['r'] = r
    U_dict['U'] = U
    U_dict['frob'] = frob
    U_dict['Uforce'] = Uforce
    U_dict['Utorque'] = Utorque

    # Lambdify all the expressions
    var_list = [Eb_sym,
                Kswt_sym,
                Krap_sym,
                rc_sym,
                Rxx, Rxy, Rxz,
                Ryx, Ryy, Ryz,
                Rzx, Rzy, Rzz,
                r_start[0],r_start[1],r_start[2],
                r_end[0],r_end[1],r_end[2],
                Sxx, Sxy, Sxz,
                Syx, Syy, Syz,
                Szx, Szy, Szz,
                Exx, Exy, Exz,
                Eyx, Eyy, Eyz,
                Ezx, Ezy, Ezz,
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
                 ori_i,ori_j,
                 Eb,
                 Kswt,Krap,
                 rc,
                 rxx,rxy,rxz,
                 ryx,ryy,ryz,
                 rzx,rzy,rzz,
                 debug=False):

    msg = "------------------------\n"
    msg += f"Computing {ti} - {tj}\n"
    if ti == "E":
        msg += "Inverting start and end\n"
        # Invert matrix R (transpose)
        rxx_c = rxx
        rxy_c = ryx
        rxz_c = rzx

        ryx_c = rxy
        ryy_c = ryy
        ryz_c = rzy

        rzx_c = rxz
        rzy_c = ryz
        rzz_c = rzz
    else:
        rxx_c = rxx
        rxy_c = rxy
        rxz_c = rxz

        ryx_c = ryx
        ryy_c = ryy
        ryz_c = ryz

        rzx_c = rzx
        rzy_c = rzy
        rzz_c = rzz

    # Compute frob numerically
    R = np.array([[rxx_c,ryx_c,rzx_c],
                  [rxy_c,ryy_c,rzy_c],
                  [rxz_c,ryz_c,rzz_c]])

    frob = np.trace(ori_j.T@ori_i@R)
    msg += f"Frobenius: {frob}\n"
    msg += f"A:\n{ori_i}\n"
    msg += f"B:\n{ori_j}\n"
    msg += f"R:\n{R}\n"
    msg += "------------------------\n"

    if debug:
        print(msg)

    analyticEval = analyticExpression(Eb,
                                      Kswt,Krap,
                                      rc,
                                      rxx_c,rxy_c,rxz_c,
                                      ryx_c,ryy_c,ryz_c,
                                      rzx_c,rzy_c,rzz_c,
                                      pos_i[0],pos_i[1],pos_i[2],
                                      pos_j[0],pos_j[1],pos_j[2],
                                      ori_i[0,0],ori_i[1,0],ori_i[2,0],
                                      ori_i[0,1],ori_i[1,1],ori_i[2,1],
                                      ori_i[0,2],ori_i[1,2],ori_i[2,2],
                                      ori_j[0,0],ori_j[1,0],ori_j[2,0],
                                      ori_j[0,1],ori_j[1,1],ori_j[2,1],
                                      ori_j[0,2],ori_j[1,2],ori_j[2,2],
                                      cn_i[0],cn_i[1],cn_i[2],
                                      cn_j[0],cn_j[1],cn_j[2])

    return analyticEval

def computePerParticleEnergyForceTorque(pos,ori,
                                        bondList,
                                        box,
                                        Eb,
                                        sigma,
                                        Kswt,Krap,
                                        rc,
                                        conn_start,conn_end,
                                        rxx,rxy,rxz,
                                        ryx,ryy,ryz,
                                        rzx,rzy,rzz):


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

    tpy2conn["S"] = conn_start
    tpy2conn["E"] = conn_end

    for i in range(N):
        p = pos[i]
        o = ori[i]

        c_S = o@conn_start + p
        c_E = o@conn_end + p

        patches.append(c_S)
        patches.append(c_E)

        patchesTpy.append("S")
        patchesTpy.append("E")

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

            if tpy_i == tpy_j:
                continue

            bondEnergy = evalAnalytic(U_dict['U'],
                                      tpy_i,tpy_j,
                                      cn_i,cn_j,
                                      pos_i,pos_j,
                                      ori_i,ori_j,
                                      Eb,
                                      Kswt, Krap,
                                      rc,
                                      rxx,rxy,rxz,
                                      ryx,ryy,ryz,
                                      rzx,rzy,rzz)

            #print(f"Energy: {bondEnergy} Frobenius: {frob_val} (ti: {tpy_i} tj: {tpy_j})")

            if bondEnergy/2.0 < minimalEnergy:
                minimalEnergy = bondEnergy/2.0
                bonds[pi] = pj
                bondsEnergy[pi] = bondEnergy/2.0


            localEnergy+= bondEnergy/2

            fs = evalAnalytic(U_dict['Uforce'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,
                              Eb,
                              Kswt, Krap,
                              rc,
                              rxx,rxy,rxz,
                              ryx,ryy,ryz,
                              rzx,rzy,rzz)

            localForce  += np.asarray([float(x) for x in fs])

            ts = evalAnalytic(U_dict['Utorque'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,
                              Eb,
                              Kswt, Krap,
                              rc,
                              rxx,rxy,rxz,
                              ryx,ryy,ryz,
                              rzx,rzy,rzz)

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





