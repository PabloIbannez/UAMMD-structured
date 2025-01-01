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

Khrm_sym  = sp.symbols("Khrm_sym")
c_i = sp.Matrix(sp.symbols('c_ix c_iy c_iz'))
c_j = sp.Matrix(sp.symbols('c_jx c_jy c_jz'))

Krap_sym  = sp.symbols("Krap_sym")
# Define the rotation matrix
Rxx, Rxy, Rxz = sp.symbols('R_xx R_xy R_xz')
Ryx, Ryy, Ryz = sp.symbols('R_yx R_yy R_yz')
Rzx, Rzy, Rzz = sp.symbols('R_zx R_zy R_zz')

R = sp.Matrix([[Rxx, Ryx, Rzx],
               [Rxy, Ryy, Rzy],
               [Rxz, Ryz, Rzz]])




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

e_x = sp.Matrix([Exx, Exy, Exz])
e_y = sp.Matrix([Eyx, Eyy, Eyz])
e_z = sp.Matrix([Ezx, Ezy, Ezz])

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

    Ub = 0.5*Khrm_sym*(r)**2 # r0 = 0

    ####################################################

    frob = (E.T@S@R).trace()
    frob = (3.0 - frob)/4.0

    Urap = Krap_sym*frob

    ####################################################

    U = Ub + Urap

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
    var_list = [Khrm_sym,
                c_i[0],c_i[1],c_i[2],
                c_j[0],c_j[1],c_j[2],
                Krap_sym,
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
                Ezx, Ezy, Ezz]

    print("Lambdifying expressions...")
    for key in U_dict:
        U_dict[key] = lambdify(var_list,U_dict[key],modules='numpy')
    print("Done.")

    return U_dict


def evalAnalytic(analyticExpression,
                 ti,tj,
                 pos_i,pos_j, # Parent
                 ori_i,ori_j,
                 Khrm,
                 cn_i,cn_j,
                 Krap,
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

    analyticEval = analyticExpression(Khrm,
                                      cn_i[0],cn_i[1],cn_i[2],
                                      cn_j[0],cn_j[1],cn_j[2],
                                      Krap,
                                      rxx_c, rxy_c, rxz_c,
                                      ryx_c, ryy_c, ryz_c,
                                      rzx_c, rzy_c, rzz_c,
                                      pos_i[0],pos_i[1],pos_i[2],
                                      pos_j[0],pos_j[1],pos_j[2],
                                      ori_i[0,0],ori_i[1,0],ori_i[2,0],
                                      ori_i[0,1],ori_i[1,1],ori_i[2,1],
                                      ori_i[0,2],ori_i[1,2],ori_i[2,2],
                                      ori_j[0,0],ori_j[1,0],ori_j[2,0],
                                      ori_j[0,1],ori_j[1,1],ori_j[2,1],
                                      ori_j[0,2],ori_j[1,2],ori_j[2,2])

    return analyticEval

def computePerParticleEnergyForceTorque(pos,ori,
                                        Khrm,
                                        conn_start,conn_end,
                                        Krap,
                                        Rxx,Rxy,Rxz,
                                        Ryx,Ryy,Ryz,
                                        Rzx,Rzy,Rzz):


    U_dict = analyticPotential()
    N = pos.shape[0]

    energies = []
    forces   = []
    torques  = []

    for i in range(N):
        energies.append(0.0)
        forces.append(np.zeros(3))
        torques.append(np.zeros(3))

    # Iterate over all consecutive pairs
    for i in range(N):
        i_prev = i - 1
        if i_prev == -1:
            i_prev = None

        i_next = i + 1
        if i_next == N:
            i_next = None

        if i_prev is not None:
            # Compute the force exerted by the previous particle
            pos_i = pos[i]
            ori_i = ori[i]

            pos_j = pos[i_prev]
            ori_j = ori[i_prev]

            conn_i = conn_end
            conn_j = conn_start

            # Compute the energy, force and torque
            energy = evalAnalytic(U_dict["U"],
                                  "E","S",
                                  pos_i,pos_j,
                                  ori_i,ori_j,
                                  Khrm,
                                  conn_i,conn_j,
                                  Krap,
                                  Rxx,Rxy,Rxz,
                                  Ryx,Ryy,Ryz,
                                  Rzx,Rzy,Rzz)

            force = evalAnalytic(U_dict["Uforce"],
                                  "E","S",
                                  pos_i,pos_j,
                                  ori_i,ori_j,
                                  Khrm,
                                  conn_i,conn_j,
                                  Krap,
                                  Rxx,Rxy,Rxz,
                                  Ryx,Ryy,Ryz,
                                  Rzx,Rzy,Rzz).flatten()

            torque = evalAnalytic(U_dict["Utorque"],
                                  "E","S",
                                  pos_i,pos_j,
                                  ori_i,ori_j,
                                  Khrm,
                                  conn_i,conn_j,
                                  Krap,
                                  Rxx,Rxy,Rxz,
                                  Ryx,Ryy,Ryz,
                                  Rzx,Rzy,Rzz).flatten()

            energies[i] += energy/2.0
            forces[i]   += force
            torques[i]  += torque

        if i_next is not None:

            pos_i = pos[i]
            ori_i = ori[i]

            pos_j = pos[i_next]
            ori_j = ori[i_next]

            conn_i = conn_start
            conn_j = conn_end

            # Compute the energy, force and torque
            energy = evalAnalytic(U_dict["U"],
                                  "S","E",
                                  pos_i,pos_j,
                                  ori_i,ori_j,
                                  Khrm,
                                  conn_i,conn_j,
                                  Krap,
                                  Rxx,Rxy,Rxz,
                                  Ryx,Ryy,Ryz,
                                  Rzx,Rzy,Rzz)

            force = evalAnalytic(U_dict["Uforce"],
                                 "S","E",
                                 pos_i,pos_j,
                                 ori_i,ori_j,
                                 Khrm,
                                 conn_i,conn_j,
                                 Krap,
                                 Rxx,Rxy,Rxz,
                                 Ryx,Ryy,Ryz,
                                 Rzx,Rzy,Rzz).flatten()



            torque = evalAnalytic(U_dict["Utorque"],
                                  "S","E",
                                  pos_i,pos_j,
                                  ori_i,ori_j,
                                  Khrm,
                                  conn_i,conn_j,
                                  Krap,
                                  Rxx,Rxy,Rxz,
                                  Ryx,Ryy,Ryz,
                                  Rzx,Rzy,Rzz).flatten()

            energies[i] += energy/2.0
            forces[i]   += force
            torques[i]  += torque

    return energies,forces,torques





