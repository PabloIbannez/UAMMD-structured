import os,sys

sys.path.append(os.path.join(sys.path[0],'..','..','scripts'))
from helixGeneration import getVx,getVy,getVz
import numpy as np
import matplotlib.pyplot as plt

import sympy as sp
# Import lambdify
from sympy.utilities.lambdify import lambdify

# Import kdtree from scipy
from MDAnalysis.lib.pkdtree import PeriodicKDTree

from collections import OrderedDict

from tqdm import tqdm

def writePairInfo(f,i,j,i_parent,j_parent,pi,pj,oi,oj,e,frc,torque):

    vxi = oi[:,0]
    vyi = oi[:,1]
    vzi = oi[:,2]

    vxj = oj[:,0]
    vyj = oj[:,1]
    vzj = oj[:,2]

    #indices_str = "id_i "+str(i)+" id_j "+str(j) + " parent_i " + str(i_parent)+" parent_j "+str(j_parent)
    indices_str = str(i)+" "+str(j) + " " + str(i_parent)+" "+str(j_parent)

    pi_str = "pos_i %.6f %.6f %.6f"%(pi[0],pi[1],pi[2])

    oi_x_str = "q_i_x %.6f %.6f %.6f"%(vxi[0],vxi[1],vxi[2])
    oi_y_str = "q_i_y %.6f %.6f %.6f"%(vyi[0],vyi[1],vyi[2])
    oi_z_str = "q_i_z %.6f %.6f %.6f"%(vzi[0],vzi[1],vzi[2])

    oi_str = oi_x_str+" "+oi_y_str+" "+oi_z_str

    pj_str = "pos_j %.6f %.6f %.6f"%(pj[0],pj[1],pj[2])

    oj_x_str = "q_j_x %.6f %.6f %.6f"%(vxj[0],vxj[1],vxj[2])
    oj_y_str = "q_j_y %.6f %.6f %.6f"%(vyj[0],vyj[1],vyj[2])
    oj_z_str = "q_j_z %.6f %.6f %.6f"%(vzj[0],vzj[1],vzj[2])

    oj_str = oj_x_str+" "+oj_y_str+" "+oj_z_str

    e_f_t_str = "%.6f %.6f %.6f %.6f %.6f %.6f %.6f"%(e,frc[0],frc[1],frc[2],torque[0],torque[1],torque[2])

    #f.write(indices_str+" "+pi_str+" "+pj_str+" "+oi_str+" "+oj_str+" "+e_f_t_str+"\n")
    f.write(indices_str+" "+e_f_t_str+"\n")

Eb_sym = sp.symbols("Eb_sym")

r_s_sym = sp.symbols("r_s_sym")
rc_sym = sp.symbols("rc_sym")

theta_start_sym = sp.symbols("theta_start_sym")
theta_end_sym   = sp.symbols("theta_end_sym")
phi_start_sym = sp.symbols("phi_start_sym")
phi_end_sym   = sp.symbols("phi_end_sym")

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
c_next = sp.Matrix(sp.symbols('c_nx c_ny c_nz'))
c_prev = sp.Matrix(sp.symbols('c_px c_py c_pz'))

# Define R matrix
#     Rxx Ryx Rzx
# R = Rxy Ryy Rzy
#     Rxz Ryz Rzz

R_sym = sp.Matrix([[sp.symbols('Rxx'),sp.symbols('Rxy'),sp.symbols('Rxz')],
                   [sp.symbols('Ryx'),sp.symbols('Ryy'),sp.symbols('Ryz')],
                   [sp.symbols('Rzx'),sp.symbols('Rzy'),sp.symbols('Rzz')]])

R_sym = R_sym.T

def analyticPotential():

    #####################################################################

    # Lab frame connections c_next_lab = S@c_next
    c_next_lab_x = s_x[0]*c_next[0] + s_y[0]*c_next[1] + s_z[0]*c_next[2]
    c_next_lab_y = s_x[1]*c_next[0] + s_y[1]*c_next[1] + s_z[1]*c_next[2]
    c_next_lab_z = s_x[2]*c_next[0] + s_y[2]*c_next[1] + s_z[2]*c_next[2]

    c_next_lab = sp.Matrix([c_next_lab_x,c_next_lab_y,c_next_lab_z])

    # Lab frame connections c_prev_lab = E@c_prev
    c_prev_lab_x = e_x[0]*c_prev[0] + e_y[0]*c_prev[1] + e_z[0]*c_prev[2]
    c_prev_lab_y = e_x[1]*c_prev[0] + e_y[1]*c_prev[1] + e_z[1]*c_prev[2]
    c_prev_lab_z = e_x[2]*c_prev[0] + e_y[2]*c_prev[1] + e_z[2]*c_prev[2]

    c_prev_lab = sp.Matrix([c_prev_lab_x,c_prev_lab_y,c_prev_lab_z])

    c_next_lab = c_next_lab + r_start
    c_prev_lab = c_prev_lab + r_end

    # Compute distance between connections

    r = sp.sqrt((c_next_lab[0]-c_prev_lab[0])**2 + (c_next_lab[1]-c_prev_lab[1])**2 + (c_next_lab[2]-c_prev_lab[2])**2)

    normalized_r = (r - r_s_sym)/(rc_sym - r_s_sym)

    Ub = sp.Piecewise((-1.0,r < r_s_sym),
                      (-sp.cos(sp.pi*normalized_r),(r >= r_s_sym) & (r <= rc_sym)),
                      ( 1.0,r > rc_sym))

    Ub = (1.0 - Ub)/2.0

    ############################################################################

    # S--R-->E
    # We have to compute S@R[:,0] which is called e_x_optimal

    cos_theta_start = sp.cos(theta_start_sym)
    cos_theta_end   = sp.cos(theta_end_sym)

    e_x_optimal = sp.Matrix([[s_x[0]*R_sym[0,0]+s_y[0]*R_sym[1,0]+s_z[0]*R_sym[2,0]],
                             [s_x[1]*R_sym[0,0]+s_y[1]*R_sym[1,0]+s_z[1]*R_sym[2,0]],
                             [s_x[2]*R_sym[0,0]+s_y[2]*R_sym[1,0]+s_z[2]*R_sym[2,0]]])

    # Now we have that cos(theta) = dot(e_x_optimal,e_x)

    cos_theta_f = e_x_optimal[0]*e_x[0] + e_x_optimal[1]*e_x[1] + e_x_optimal[2]*e_x[2]

    normalized_theta_f = (cos_theta_f - cos_theta_start)/(cos_theta_end - cos_theta_start)

    Ua_forward = sp.Piecewise((-1.0, cos_theta_f > cos_theta_start),
                              (-sp.cos(sp.pi * normalized_theta_f), (cos_theta_f <= cos_theta_start) & (cos_theta_f >= cos_theta_end)),
                              (1.0, cos_theta_f < cos_theta_end))

    Ua_forward = (1.0-Ua_forward)/2.0

    # Now we compute Theta_B potential
    # S<--R.T--E
    # We have to compute E@R.T[:,0] which is called s_x_optimal

    s_x_optimal = sp.Matrix([[e_x[0]*R_sym[0,0]+e_y[0]*R_sym[0,1]+e_z[0]*R_sym[0,2]],
                             [e_x[1]*R_sym[0,0]+e_y[1]*R_sym[0,1]+e_z[1]*R_sym[0,2]],
                             [e_x[2]*R_sym[0,0]+e_y[2]*R_sym[0,1]+e_z[2]*R_sym[0,2]]])

    ## Now we have that cos(theta) = dot(s_x_optimal,s_x)

    cos_theta_b = s_x_optimal[0]*s_x[0] + s_x_optimal[1]*s_x[1] + s_x_optimal[2]*s_x[2]

    normalized_theta_b = (cos_theta_b - cos_theta_start)/(cos_theta_end - cos_theta_start)

    Ua_backward = sp.Piecewise((-1.0,cos_theta_b>cos_theta_start),
                               (-sp.cos(sp.pi*normalized_theta_b),(cos_theta_b<=cos_theta_start) & (cos_theta_b>=cos_theta_end)),
                               ( 1.0,cos_theta_b<cos_theta_end))

    Ua_backward = (1.0-Ua_backward)/2.0

    Ua = (Ua_forward+Ua_backward)/2.0

    #####################################################

    # Now we comute forward and backward potential for phi

    cos_phi_start = sp.cos(phi_start_sym)
    cos_phi_end   = sp.cos(phi_end_sym)

    e_y_optimal = sp.Matrix([[s_x[0]*R_sym[0,1]+s_y[0]*R_sym[1,1]+s_z[0]*R_sym[2,1]],
                             [s_x[1]*R_sym[0,1]+s_y[1]*R_sym[1,1]+s_z[1]*R_sym[2,1]],
                             [s_x[2]*R_sym[0,1]+s_y[2]*R_sym[1,1]+s_z[2]*R_sym[2,1]]])

    cos_phi_f = e_y_optimal[0]*e_y[0] + e_y_optimal[1]*e_y[1] + e_y_optimal[2]*e_y[2]

    normalized_phi_f = (cos_phi_f - cos_phi_start)/(cos_phi_end - cos_phi_start)

    Ud_forward = sp.Piecewise((-1.0,cos_phi_f>cos_phi_start),
                              (-sp.cos(sp.pi*normalized_phi_f),(cos_phi_f<=cos_phi_start) & (cos_phi_f>=cos_phi_end)),
                              ( 1.0,cos_phi_f<cos_phi_end))

    Ud_forward = (1.0-Ud_forward)/2.0

    # Backward potential

    s_y_optimal = sp.Matrix([[e_x[0]*R_sym[1,0]+e_y[0]*R_sym[1,1]+e_z[0]*R_sym[1,2]],
                             [e_x[1]*R_sym[1,0]+e_y[1]*R_sym[1,1]+e_z[1]*R_sym[1,2]],
                             [e_x[2]*R_sym[1,0]+e_y[2]*R_sym[1,1]+e_z[2]*R_sym[1,2]]])

    cos_phi_b = s_y_optimal[0]*s_y[0] + s_y_optimal[1]*s_y[1] + s_y_optimal[2]*s_y[2]

    normalized_phi_b = (cos_phi_b - cos_phi_start)/(cos_phi_end - cos_phi_start)

    Ud_backward = sp.Piecewise((-1.0,cos_phi_b>cos_phi_start),
                               (-sp.cos(sp.pi*normalized_phi_b),(cos_phi_b<=cos_phi_start) & (cos_phi_b>=cos_phi_end)),
                               ( 1.0,cos_phi_b<cos_phi_end))

    Ud_backward = (1.0-Ud_backward)/2.0

    Ud = (Ud_forward+Ud_backward)/2.0

    ####################################################

    U = -Eb_sym*Ub*Ua*Ud

    #print("Computing S force ...")
    UforceS   = -U.diff(r_start)
    #print("Computing E force ...")
    UforceE   = -U.diff(r_end)
    #print("Computing S torque ...")
    UtorqueS  = - s_x.cross(U.diff(s_x)) - s_y.cross(U.diff(s_y)) - s_z.cross(U.diff(s_z))
    #print("Computing E torque ...")
    UtorqueE  = - e_x.cross(U.diff(e_x)) - e_y.cross(U.diff(e_y)) - e_z.cross(U.diff(e_z))

    #dUads_x    = Ua.diff(s_x)
    #dUads_y    = Ua.diff(s_y)
    #dUads_z    = Ua.diff(s_z)

    #dUafds_x   = Ua_forward.diff(s_x)
    #dUafds_y   = Ua_forward.diff(s_y)
    #dUafds_z   = Ua_forward.diff(s_z)

    #dUabds_x   = Ua_backward.diff(s_x)
    #dUabds_y   = Ua_backward.diff(s_y)
    #dUabds_z   = Ua_backward.diff(s_z)

    #dUade_x    = Ua.diff(e_x)
    #dUade_y    = Ua.diff(e_y)
    #dUade_z    = Ua.diff(e_z)

    #dUafde_x   = Ua_forward.diff(e_x)
    #dUafde_y   = Ua_forward.diff(e_y)
    #dUafde_z   = Ua_forward.diff(e_z)

    #dUabde_x   = Ua_backward.diff(e_x)
    #dUabde_y   = Ua_backward.diff(e_y)
    #dUabde_z   = Ua_backward.diff(e_z)

    ###################################

    #dUdds_x    = Ud.diff(s_x)
    #dUdds_y    = Ud.diff(s_y)
    #dUdds_z    = Ud.diff(s_z)

    #dUdfds_x    = Ud_forward.diff(s_x)
    #dUdfds_y    = Ud_forward.diff(s_y)
    #dUdfds_z    = Ud_forward.diff(s_z)

    #dUdbds_x    = Ud_backward.diff(s_x)
    #dUdbds_y    = Ud_backward.diff(s_y)
    #dUdbds_z    = Ud_backward.diff(s_z)

    #dUdde_x    = Ud.diff(e_x)
    #dUdde_y    = Ud.diff(e_y)
    #dUdde_z    = Ud.diff(e_z)

    #dUdfde_x    = Ud_forward.diff(e_x)
    #dUdfde_y    = Ud_forward.diff(e_y)
    #dUdfde_z    = Ud_forward.diff(e_z)

    #dUdbde_x    = Ud_backward.diff(e_x)
    #dUdbde_y    = Ud_backward.diff(e_y)
    #dUdbde_z    = Ud_backward.diff(e_z)

    ##################################

    U_dict = OrderedDict()

    U_dict['r'] = r
    U_dict['U'] = U
    U_dict['UforceS'] = UforceS
    U_dict['UforceE'] = UforceE
    U_dict['UtorqueS'] = UtorqueS
    U_dict['UtorqueE'] = UtorqueE

    #U_dict['Ub'] = -Eb_sym*Ub
    #U_dict['Uo'] =  Ua*Ud

    #U_dict['Ua'] = Ua

    #U_dict['Ua_forward']  = Ua_forward
    #U_dict['Ua_backward'] = Ua_backward

    #U_dict['Ud'] = Ud

    #U_dict['Ud_forward']  = Ud_forward
    #U_dict['Ud_backward'] = Ud_backward

    #U_dict['dUads_x']  = dUads_x
    #U_dict['dUads_y']  = dUads_y
    #U_dict['dUads_z']  = dUads_z

    #U_dict['dUafds_x'] = dUafds_x
    #U_dict['dUafds_y'] = dUafds_y
    #U_dict['dUafds_z'] = dUafds_z

    #U_dict['dUabds_x'] = dUabds_x
    #U_dict['dUabds_y'] = dUabds_y
    #U_dict['dUabds_z'] = dUabds_z

    #U_dict['dUdds_x']  = dUdds_x
    #U_dict['dUdds_y']  = dUdds_y
    #U_dict['dUdds_z']  = dUdds_z

    #U_dict['dUdfds_x'] = dUdfds_x
    #U_dict['dUdfds_y'] = dUdfds_y
    #U_dict['dUdfds_z'] = dUdfds_z

    #U_dict['dUdbds_x'] = dUdbds_x
    #U_dict['dUdbds_y'] = dUdbds_y
    #U_dict['dUdbds_z'] = dUdbds_z

    #U_dict['dUade_x']  = dUade_x
    #U_dict['dUade_y']  = dUade_y
    #U_dict['dUade_z']  = dUade_z

    #U_dict['dUafde_x'] = dUafde_x
    #U_dict['dUafde_y'] = dUafde_y
    #U_dict['dUafde_z'] = dUafde_z

    #U_dict['dUabde_x'] = dUabde_x
    #U_dict['dUabde_y'] = dUabde_y
    #U_dict['dUabde_z'] = dUabde_z

    #U_dict['dUdde_x']  = dUdde_x
    #U_dict['dUdde_y']  = dUdde_y
    #U_dict['dUdde_z']  = dUdde_z

    #U_dict['dUdfde_x'] = dUdfde_x
    #U_dict['dUdfde_y'] = dUdfde_y
    #U_dict['dUdfde_z'] = dUdfde_z

    #U_dict['dUdbde_x'] = dUdbde_x
    #U_dict['dUdbde_y'] = dUdbde_y
    #U_dict['dUdbde_z'] = dUdbde_z

    # Lambdify all the expressions
    var_list = [Eb_sym,
                r_s_sym,
                rc_sym,
                theta_start_sym,theta_end_sym,
                phi_start_sym,phi_end_sym,
                r_start[0],r_start[1],r_start[2],
                r_end[0],r_end[1],r_end[2],
                s_x[0],s_x[1],s_x[2],
                s_y[0],s_y[1],s_y[2],
                s_z[0],s_z[1],s_z[2],
                e_x[0],e_x[1],e_x[2],
                e_y[0],e_y[1],e_y[2],
                e_z[0],e_z[1],e_z[2],
                c_next[0],c_next[1],c_next[2],
                c_prev[0],c_prev[1],c_prev[2],
                R_sym[0,0],R_sym[0,1],R_sym[0,2],
                R_sym[1,0],R_sym[1,1],R_sym[1,2],
                R_sym[2,0],R_sym[2,1],R_sym[2,2]]

    print("Lambdifying expressions...")
    for key in U_dict:
        U_dict[key] = lambdify(var_list,U_dict[key],modules='numpy')
    print("Done.")

    return U_dict


def evalAnalytic(analyticExpression,
                 ti,tj,
                 cn_i,cn_j,
                 pos_i,pos_j, # Parent
                 ori_i,ori_j,R_H,Eb,
                 r_s,rc,
                 theta_start,theta_end,
                 phi_start,phi_end):

    cn_s  = cn_i
    pos_s = pos_i
    ori_s = ori_i

    cn_e  = cn_j
    pos_e = pos_j
    ori_e = ori_j

    if ti == "S" and tj == "E":
        R = R_H
    elif ti == "E" and tj == "S":
        R = R_H.T
    else:
        print("ti: ",ti," tj: ",tj)
        raise ValueError("Wrong type of connection")

    analyticEval = analyticExpression(Eb,
                                      r_s,rc,
                                      theta_start,theta_end,
                                      phi_start,phi_end,
                                      pos_s[0],pos_s[1],pos_s[2],
                                      pos_e[0],pos_e[1],pos_e[2],
                                      ori_s[0,0],ori_s[1,0],ori_s[2,0],
                                      ori_s[0,1],ori_s[1,1],ori_s[2,1],
                                      ori_s[0,2],ori_s[1,2],ori_s[2,2],
                                      ori_e[0,0],ori_e[1,0],ori_e[2,0],
                                      ori_e[0,1],ori_e[1,1],ori_e[2,1],
                                      ori_e[0,2],ori_e[1,2],ori_e[2,2],
                                      cn_s[0],cn_s[1],cn_s[2],
                                      cn_e[0],cn_e[1],cn_e[2],
                                      R[0,0],R[0,1],R[0,2],
                                      R[1,0],R[1,1],R[1,2],
                                      R[2,0],R[2,1],R[2,2])

    return analyticEval

def computePerParticleHelixEnergyForceTorque(pos,ori,
                                             bondList,
                                             box,
                                             Eb,
                                             sigma,
                                             r_s,rc,
                                             theta_start,theta_end,
                                             phi_start,phi_end,
                                             R_H,conn_next,conn_prev):

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
    patchesTpy = []
    parentId   = []

    #test = open("test.sp","w")
    for i in range(N):
        p = pos[i]
        o = ori[i]

        c_next_pos = o@conn_next + p
        c_prev_pos = o@conn_prev + p

        patches.append(c_prev_pos)
        patches.append(c_next_pos)

        patchesTpy.append("E")
        patchesTpy.append("S")

        parentId.append(i)
        parentId.append(i)

        #format x y z sigma/20 c
        #test.write(str(p[0])+" "+str(p[1])+" "+str(p[2])+" "+str(sigma/2)+" 0\n")
        #test.write(str(c_next_pos[0])+" "+str(c_next_pos[1])+" "+str(c_next_pos[2])+" "+str(sigma/10)+" 1\n")
        #test.write(str(c_prev_pos[0])+" "+str(c_prev_pos[1])+" "+str(c_prev_pos[2])+" "+str(sigma/10)+" 2\n")

    patch2parent = {}
    for i in range(len(patches)):
        patch2parent[i] = parentId[i]

    # Create KDtree using patches and rc

    tree = PeriodicKDTree(box=np.asarray(box+[90,90,90],dtype=np.float32))
    tree.set_coords(np.asarray(patches,dtype=np.float32),cutoff=rc)

    # Query pairs
    pairs = tree.search_pairs(rc)

    # Remove all pairs which types are the same
    pairs = [p for p in pairs if patchesTpy[p[0]] != patchesTpy[p[1]]]

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

        if tpy_i == "S":
            cn_i = conn_next
        elif tpy_i == "E":
            cn_i = conn_prev
        else:
            raise ValueError("Wrong type of patch")

        localEnergy = 0.0

        localForce  = np.zeros(3)
        localTorque = np.zeros(3)

        minimalEnergy = 0.0

        for pj in neigh:

            pos_j = pos[patch2parent[pj]]
            ori_j = ori[patch2parent[pj]]

            tpy_j = patchesTpy[pj]

            if tpy_j == "S":
                cn_j = conn_next
            elif tpy_j == "E":
                cn_j = conn_prev
            else:
                raise ValueError("Wrong type of patch")

            if tpy_i == tpy_j:
                raise ValueError("Error: Same type of patch")

            bondEnergy = evalAnalytic(U_dict['U'],
                                      tpy_i,tpy_j,
                                      cn_i,cn_j,
                                      pos_i,pos_j,
                                      ori_i,ori_j,R_H,Eb,
                                      r_s,rc,
                                      theta_start,theta_end,
                                      phi_start,phi_end)

            if bondEnergy/2.0 < minimalEnergy:
                minimalEnergy = bondEnergy/2.0
                bonds[pi] = pj
                bondsEnergy[pi] = bondEnergy/2.0

            #r = evalAnalytic(U_dict['r'],
            #                 tpy_i,tpy_j,
            #                 cn_i,cn_j,
            #                 pos_i,pos_j,
            #                 ori_i,ori_j,R_H,Eb,rc,Kb,Ka,Kd)

            #print(patch2parent[pi],patch2parent[pj],bondEnergy,r)

            localEnergy+= bondEnergy/2

            fs = evalAnalytic(U_dict['UforceS'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,R_H,Eb,
                              r_s,rc,
                              theta_start,theta_end,
                              phi_start,phi_end)

            localForce  += np.asarray([float(x) for x in fs])

            ts = evalAnalytic(U_dict['UtorqueS'],
                              tpy_i,tpy_j,
                              cn_i,cn_j,
                              pos_i,pos_j,
                              ori_i,ori_j,R_H,Eb,
                              r_s,rc,
                              theta_start,theta_end,
                              phi_start,phi_end)

            localTorque += np.asarray([float(x) for x in ts])

            #writePairInfo(logFile,pi,pj,patch2parent[pi],patch2parent[pj],patches[pi],patches[pj],ori[patch2parent[pi]],ori[patch2parent[pj]],bondEnergy/2.0,fs,ts)

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





