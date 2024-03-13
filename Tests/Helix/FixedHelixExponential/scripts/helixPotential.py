import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.utilities.lambdify import lambdify

from collections import OrderedDict

from tqdm import tqdm

Eb_sym = sp.symbols("Eb_sym")

Kb_sym = sp.symbols("Kb_sym")
Ka_sym = sp.symbols("Ka_sym")
Kd_sym = sp.symbols("Kd_sym")

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

    # The potential is 0.5*Kb*r^2

    Ub = 0.5*Kb_sym*r**2

    # Theta_F potential is (exp(-Ka(1.0-cos(theta)))-A)/(B-A) where A = exp(-2Ka) and B = 1.0
    # S--R-->E
    # We have to compute S@R[:,0] which is called e_x_optimal

    e_x_optimal = sp.Matrix([[s_x[0]*R_sym[0,0]+s_y[0]*R_sym[1,0]+s_z[0]*R_sym[2,0]],
                             [s_x[1]*R_sym[0,0]+s_y[1]*R_sym[1,0]+s_z[1]*R_sym[2,0]],
                             [s_x[2]*R_sym[0,0]+s_y[2]*R_sym[1,0]+s_z[2]*R_sym[2,0]]])

    # Now we have that cos(theta) = dot(e_x_optimal,e_x)

    cos_theta_f = e_x_optimal[0]*e_x[0] + e_x_optimal[1]*e_x[1] + e_x_optimal[2]*e_x[2]

    Ua_forward = (sp.exp(-Ka_sym*(1.0-cos_theta_f))-sp.exp(-2.0*Ka_sym))/(1.0-sp.exp(-2.0*Ka_sym))

    # Now we compute Theta_B potential
    # S<--R.T--E
    # We have to compute E@R.T[:,0] which is called s_x_optimal

    s_x_optimal = sp.Matrix([[e_x[0]*R_sym[0,0]+e_y[0]*R_sym[0,1]+e_z[0]*R_sym[0,2]],
                             [e_x[1]*R_sym[0,0]+e_y[1]*R_sym[0,1]+e_z[1]*R_sym[0,2]],
                             [e_x[2]*R_sym[0,0]+e_y[2]*R_sym[0,1]+e_z[2]*R_sym[0,2]]])

    ## Now we have that cos(theta) = dot(s_x_optimal,s_x)

    cos_theta_b = s_x_optimal[0]*s_x[0] + s_x_optimal[1]*s_x[1] + s_x_optimal[2]*s_x[2]

    Ua_backward = (sp.exp(-Ka_sym*(1.0-cos_theta_b))-sp.exp(-2.0*Ka_sym))/(1.0-sp.exp(-2.0*Ka_sym))

    Ua = (Ua_forward+Ua_backward)/2.0

    #####################################################

    # Now we comute forward and backward potential for phi

    e_y_optimal = sp.Matrix([[s_x[0]*R_sym[0,1]+s_y[0]*R_sym[1,1]+s_z[0]*R_sym[2,1]],
                             [s_x[1]*R_sym[0,1]+s_y[1]*R_sym[1,1]+s_z[1]*R_sym[2,1]],
                             [s_x[2]*R_sym[0,1]+s_y[2]*R_sym[1,1]+s_z[2]*R_sym[2,1]]])

    cos_phi_f = e_y_optimal[0]*e_y[0] + e_y_optimal[1]*e_y[1] + e_y_optimal[2]*e_y[2]

    Ud_forward = (sp.exp(-Kd_sym*(1.0-cos_phi_f))-sp.exp(-2.0*Kd_sym))/(1.0-sp.exp(-2.0*Kd_sym))

    # Backward potential

    s_y_optimal = sp.Matrix([[e_x[0]*R_sym[1,0]+e_y[0]*R_sym[1,1]+e_z[0]*R_sym[1,2]],
                             [e_x[1]*R_sym[1,0]+e_y[1]*R_sym[1,1]+e_z[1]*R_sym[1,2]],
                             [e_x[2]*R_sym[1,0]+e_y[2]*R_sym[1,1]+e_z[2]*R_sym[1,2]]])

    cos_phi_b = s_y_optimal[0]*s_y[0] + s_y_optimal[1]*s_y[1] + s_y_optimal[2]*s_y[2]

    Ud_backward = (sp.exp(-Kd_sym*(1.0-cos_phi_b))-sp.exp(-2.0*Kd_sym))/(1.0-sp.exp(-2.0*Kd_sym))

    Ud = (Ud_forward+Ud_backward)/2.0

    ####################################################

    U = -Eb_sym*Ua*Ud + Ub

    UforceS   = -U.diff(r_start)
    UforceE   = -U.diff(r_end)

    dUads_x    = Ua.diff(s_x)
    dUads_y    = Ua.diff(s_y)
    dUads_z    = Ua.diff(s_z)

    dUafds_x   = Ua_forward.diff(s_x)
    dUafds_y   = Ua_forward.diff(s_y)
    dUafds_z   = Ua_forward.diff(s_z)

    dUabds_x   = Ua_backward.diff(s_x)
    dUabds_y   = Ua_backward.diff(s_y)
    dUabds_z   = Ua_backward.diff(s_z)

    dUade_x    = Ua.diff(e_x)
    dUade_y    = Ua.diff(e_y)
    dUade_z    = Ua.diff(e_z)

    dUafde_x   = Ua_forward.diff(e_x)
    dUafde_y   = Ua_forward.diff(e_y)
    dUafde_z   = Ua_forward.diff(e_z)

    dUabde_x   = Ua_backward.diff(e_x)
    dUabde_y   = Ua_backward.diff(e_y)
    dUabde_z   = Ua_backward.diff(e_z)

    ##################################

    dUdds_x    = Ud.diff(s_x)
    dUdds_y    = Ud.diff(s_y)
    dUdds_z    = Ud.diff(s_z)

    dUdfds_x    = Ud_forward.diff(s_x)
    dUdfds_y    = Ud_forward.diff(s_y)
    dUdfds_z    = Ud_forward.diff(s_z)

    dUdbds_x    = Ud_backward.diff(s_x)
    dUdbds_y    = Ud_backward.diff(s_y)
    dUdbds_z    = Ud_backward.diff(s_z)

    dUdde_x    = Ud.diff(e_x)
    dUdde_y    = Ud.diff(e_y)
    dUdde_z    = Ud.diff(e_z)

    dUdfde_x    = Ud_forward.diff(e_x)
    dUdfde_y    = Ud_forward.diff(e_y)
    dUdfde_z    = Ud_forward.diff(e_z)

    dUdbde_x    = Ud_backward.diff(e_x)
    dUdbde_y    = Ud_backward.diff(e_y)
    dUdbde_z    = Ud_backward.diff(e_z)

    ##################################

    UtorqueS  = - s_x.cross(U.diff(s_x)) - s_y.cross(U.diff(s_y)) - s_z.cross(U.diff(s_z))
    UtorqueE  = - e_x.cross(U.diff(e_x)) - e_y.cross(U.diff(e_y)) - e_z.cross(U.diff(e_z))

    U_dict = OrderedDict()

    U_dict['U'] = U
    U_dict['UforceS'] = UforceS
    U_dict['UforceE'] = UforceE

    U_dict['Ua'] = Ua

    U_dict['Ua_forward']  = Ua_forward
    U_dict['Ua_backward'] = Ua_backward

    U_dict['Ud'] = Ud

    U_dict['Ud_forward']  = Ud_forward
    U_dict['Ud_backward'] = Ud_backward

    U_dict['dUads_x']  = dUads_x
    U_dict['dUads_y']  = dUads_y
    U_dict['dUads_z']  = dUads_z

    U_dict['dUafds_x'] = dUafds_x
    U_dict['dUafds_y'] = dUafds_y
    U_dict['dUafds_z'] = dUafds_z

    U_dict['dUabds_x'] = dUabds_x
    U_dict['dUabds_y'] = dUabds_y
    U_dict['dUabds_z'] = dUabds_z

    U_dict['dUdds_x']  = dUdds_x
    U_dict['dUdds_y']  = dUdds_y
    U_dict['dUdds_z']  = dUdds_z

    U_dict['dUdfds_x'] = dUdfds_x
    U_dict['dUdfds_y'] = dUdfds_y
    U_dict['dUdfds_z'] = dUdfds_z

    U_dict['dUdbds_x'] = dUdbds_x
    U_dict['dUdbds_y'] = dUdbds_y
    U_dict['dUdbds_z'] = dUdbds_z

    U_dict['dUade_x']  = dUade_x
    U_dict['dUade_y']  = dUade_y
    U_dict['dUade_z']  = dUade_z

    U_dict['dUafde_x'] = dUafde_x
    U_dict['dUafde_y'] = dUafde_y
    U_dict['dUafde_z'] = dUafde_z

    U_dict['dUabde_x'] = dUabde_x
    U_dict['dUabde_y'] = dUabde_y
    U_dict['dUabde_z'] = dUabde_z

    U_dict['dUdde_x']  = dUdde_x
    U_dict['dUdde_y']  = dUdde_y
    U_dict['dUdde_z']  = dUdde_z

    U_dict['dUdfde_x'] = dUdfde_x
    U_dict['dUdfde_y'] = dUdfde_y
    U_dict['dUdfde_z'] = dUdfde_z

    U_dict['dUdbde_x'] = dUdbde_x
    U_dict['dUdbde_y'] = dUdbde_y
    U_dict['dUdbde_z'] = dUdbde_z

    U_dict['UtorqueS'] = UtorqueS
    U_dict['UtorqueE'] = UtorqueE

    # Lambdify all the expressions
    var_list = [Eb_sym,
                Kb_sym,Ka_sym,Kd_sym,
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
                 ori_i,ori_j,R_H,Eb,Kb,Ka,Kd):

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

    #subs_dict = {Eb_sym:Eb,
    #             Kb_sym:Kb,Ka_sym:Ka,Kd_sym:Kd,
    #             r_start[0]:pos_s[0],r_start[1]:pos_s[1],r_start[2]:pos_s[2],
    #             r_end[0]:pos_e[0],r_end[1]:pos_e[1],r_end[2]:pos_e[2],
    #             s_x[0]:ori_s[0,0],s_x[1]:ori_s[1,0],s_x[2]:ori_s[2,0],
    #             s_y[0]:ori_s[0,1],s_y[1]:ori_s[1,1],s_y[2]:ori_s[2,1],
    #             s_z[0]:ori_s[0,2],s_z[1]:ori_s[1,2],s_z[2]:ori_s[2,2],
    #             e_x[0]:ori_e[0,0],e_x[1]:ori_e[1,0],e_x[2]:ori_e[2,0],
    #             e_y[0]:ori_e[0,1],e_y[1]:ori_e[1,1],e_y[2]:ori_e[2,1],
    #             e_z[0]:ori_e[0,2],e_z[1]:ori_e[1,2],e_z[2]:ori_e[2,2],
    #             c_next[0]:cn_s[0],c_next[1]:cn_s[1],c_next[2]:cn_s[2],
    #             c_prev[0]:cn_e[0],c_prev[1]:cn_e[1],c_prev[2]:cn_e[2],
    #             R_sym[0,0]:R[0,0],R_sym[0,1]:R[0,1],R_sym[0,2]:R[0,2],
    #             R_sym[1,0]:R[1,0],R_sym[1,1]:R[1,1],R_sym[1,2]:R[1,2],
    #             R_sym[2,0]:R[2,0],R_sym[2,1]:R[2,1],R_sym[2,2]:R[2,2]}

    #analyticEval = analyticExpression.evalf(20,subs=subs_dict)

    analyticEval = analyticExpression(Eb,
                                      Kb,Ka,Kd,
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

def computePerParticleHelixEnergyForceTorque(pos,ori,Eb,Kb,Ka,Kd,R_H,conn_next,conn_prev):

    U_dict = analyticPotential()
    N = pos.shape[0]

    energy  = []
    forces  = []
    torques = []

    for i in range(N):
        energy.append(0.0)
        forces.append(np.zeros(3))
        torques.append(np.zeros(3))

    for i in tqdm(range(N-1)):
        pos_i = pos[i]
        pos_j = pos[i+1]

        ori_i = ori[i]
        ori_j = ori[i+1]

        cn_i = conn_next
        cn_j = conn_prev

        bondEnergy = evalAnalytic(U_dict['U'],
                                  "S","E",
                                  cn_i,cn_j,
                                  pos_i,pos_j,
                                  ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        energy[i]   += bondEnergy/2.0
        energy[i+1] += bondEnergy/2.0

        fs = evalAnalytic(U_dict['UforceS'],
                          "S","E",
                          cn_i,cn_j,
                          pos_i,pos_j,
                          ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        fe = evalAnalytic(U_dict['UforceE'],
                          "S","E",
                          cn_i,cn_j,
                          pos_i,pos_j,
                          ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        forces[i]   += np.asarray([float(x) for x in fs])
        forces[i+1] += np.asarray([float(x) for x in fe])

        ts = evalAnalytic(U_dict['UtorqueS'],
                          "S","E",
                          cn_i,cn_j,
                          pos_i,pos_j,
                          ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        te = evalAnalytic(U_dict['UtorqueE'],
                          "S","E",
                          cn_i,cn_j,
                          pos_i,pos_j,
                          ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        torques[i]   += np.asarray([float(x) for x in ts])
        torques[i+1] += np.asarray([float(x) for x in te])

        #for entry,expression in U_dict.items():
        #    expression_v = evalAnalytic(expression,
        #                                "S","E",
        #                                cn_i,cn_j,
        #                                pos_i,pos_j,
        #                                ori_i,ori_j,R_H,Eb,Kb,Ka,Kd)

        #    if "ds_" not in entry and "de_" not in entry:
        #        index = ""
        #    else:
        #        if "de_" in entry:
        #            index = i + 1
        #        else:
        #            index = i

        #    try:
        #        print(i,index,entry,list(expression_v))
        #    except:
        #        print(i,index,entry,expression_v)

        #print("####################################")

    return energy,forces,torques





