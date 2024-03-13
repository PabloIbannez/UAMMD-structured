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

E_sym = sp.symbols("E_sym")
D_sym = sp.symbols("D_sym")

# Define position vector for the paricle
r_linker_sym = sp.Matrix(sp.symbols('r_lx r_ly r_lz'))

surf_pos_sym = sp.symbols("surf_pos")

# Define orientation vectors for the particle
o_x = sp.Matrix(sp.symbols('o_xx o_xy o_xz'))
o_y = sp.Matrix(sp.symbols('o_yx o_yy o_yz'))
o_z = sp.Matrix(sp.symbols('o_zx o_zy o_zz'))

linker_pos_sym = sp.Matrix(sp.symbols('lp_x lp_y lp_z'))

def analyticPotential():

    #####################################################################

    # Lab frame for linker_lab = O@linker_pos
    linker_lab_x = o_x[0]*linker_pos_sym[0] + o_y[0]*linker_pos_sym[1] + o_z[0]*linker_pos_sym[2]
    linker_lab_y = o_x[1]*linker_pos_sym[0] + o_y[1]*linker_pos_sym[1] + o_z[1]*linker_pos_sym[2]
    linker_lab_z = o_x[2]*linker_pos_sym[0] + o_y[2]*linker_pos_sym[1] + o_z[2]*linker_pos_sym[2]

    linker_lab = sp.Matrix([linker_lab_x,linker_lab_y,linker_lab_z])

    linker_lab = linker_lab + r_linker_sym

    # Compute distance between linker and surface

    dz = linker_lab[2] - surf_pos_sym

    # Compute the potential U = -E*exp(-dz*dz/(2.0*D))

    U = -E_sym*sp.exp(-dz*dz/(2.0*D_sym))

    ####################################################

    Uforce   = -U.diff(r_linker_sym)
    Utorque  = - o_x.cross(U.diff(o_x)) - o_y.cross(U.diff(o_y)) - o_z.cross(U.diff(o_z))

    U_dict = OrderedDict()

    U_dict['dz'] = dz
    U_dict['r']  = linker_lab

    U_dict['U'] = U
    U_dict['Uforce'] = Uforce
    U_dict['Utorque'] = Utorque

    # Lambdify all the expressions
    var_list = [E_sym,
                D_sym,
                r_linker_sym[0],r_linker_sym[1],r_linker_sym[2],
                o_x[0],o_x[1],o_x[2],
                o_y[0],o_y[1],o_y[2],
                o_z[0],o_z[1],o_z[2],
                linker_pos_sym[0],linker_pos_sym[1],linker_pos_sym[2],
                surf_pos_sym]

    print("Lambdifying expressions...")
    for key in U_dict:
        U_dict[key] = lambdify(var_list,U_dict[key],modules='numpy')
    print("Done.")

    return U_dict


def evalAnalytic(analyticExpression,
                 linker_pos,
                 pos, # Parent
                 ori,
                 E,D,surf_pos):

    #print("Evaluting expression ...")

    #analyticEval = analyticExpression.evalf(20,subs=subs_dict)

    analyticEval = analyticExpression(E,D,
                                      pos[0],pos[1],pos[2],
                                      ori[0][0],ori[1][0],ori[2][0],
                                      ori[0][1],ori[1][1],ori[2][1],
                                      ori[0][2],ori[1][2],ori[2][2],
                                      linker_pos[0],linker_pos[1],linker_pos[2],
                                      surf_pos)

    return analyticEval

def computePerParticleLinkerEnergyForceTorque(pos,ori,
                                              E,D,linker_pos,surf_pos):

    #print("R_H: ",R_H)

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
    parentId   = []

    #test = open("test.sp","w")
    for i in range(N):
        p = pos[i]
        o = ori[i]

        l = o@linker_pos + p

        patches.append(l)
        parentId.append(i)

    patch2parent = {}
    for i in range(len(patches)):
        patch2parent[i] = parentId[i]

    for i,p in enumerate(patches):

        pos_i = pos[patch2parent[i]]
        ori_i = ori[patch2parent[i]]

        energy[i] = evalAnalytic(U_dict['U'],
                                 linker_pos,
                                 pos_i, # Parent
                                 ori_i,
                                 E,D,surf_pos)

        fs = evalAnalytic(U_dict['Uforce'],
                          linker_pos,
                          pos_i, # Parent
                          ori_i,
                          E,D,surf_pos)

        forces[i] = np.asarray([float(x) for x in fs])

        ts = evalAnalytic(U_dict['Utorque'],
                          linker_pos,
                          pos_i, # Parent
                          ori_i,
                          E,D,surf_pos)

        torques[i] = np.asarray([float(x) for x in ts])

        dz = evalAnalytic(U_dict['dz'],
                          linker_pos,
                          pos_i, # Parent
                          ori_i,
                          E,D,surf_pos)

        r = evalAnalytic(U_dict['r'],
                         linker_pos,
                         pos_i, # Parent
                         ori_i,
                         E,D,surf_pos)

        print(r,dz)


    return energy,forces,torques





