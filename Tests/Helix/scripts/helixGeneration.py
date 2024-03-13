import numpy as np

from scipy.optimize import fsolve
from scipy.spatial.transform import Rotation

import itertools

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

def smallRandomRotation(angleRange):

    # Generate a random unit vector
    random_vector = np.random.rand(3)
    random_vector /= np.linalg.norm(random_vector)

    # Generate a random angle between -angleRange and angleRange
    random_angle = np.random.rand()*angleRange*2.0-angleRange

    # Create the rotation vector (axis * angle)
    rotation_vector = random_vector * random_angle

    # Create the rotation
    rotation = Rotation.from_rotvec(rotation_vector)

    return rotation

def helixEquation(s,a,b,e):

    H = np.sqrt(a**2 + b**2)

    x = a*np.cos(s/H)
    y = a*np.sin(s/H)*e
    z = b*s/H

    return np.asarray([x,y,z])

def helixPointsDistanceMinimization(s1,s2,a,b,e,sigma):

    p1 = helixEquation(s1,a,b,e)
    p2 = helixEquation(s2,a,b,e)

    dx = p1[0]-p2[0]
    dy = p1[1]-p2[1]
    dz = p1[2]-p2[2]

    d = np.sqrt(dx*dx+dy*dy+dz*dz)-sigma

    return d

def computeHelixParticlesDistance(a,b,e,sigma):

    try:
        sOne = fsolve(helixPointsDistanceMinimization,sigma,args=(0.0,a,b,e,sigma),xtol=1e-5)[0] # x0 = sigma
        #Check if sOne is correct
        diff = np.abs(helixPointsDistanceMinimization(sOne,0.0,a,b,e,sigma))
        if diff >1e-5:
            print(f"Error computing sOne, diff: {diff}")
            raise Exception("sOne is not correct")
    except:
        print(f"Error computing sOne")
        raise Exception("Error computing sOne")

    return sOne


def discreteHelixFrame(i,sOne,a,b,e):

    pos     = helixEquation( i*sOne,a,b,e)
    posNext = helixEquation((i+1)*sOne,a,b,e)

    dr = posNext-pos
    ex = dr/np.linalg.norm(dr)

    ey = np.cross(np.asarray([0.0,0.0,1.0]),ex)
    ey = ey/np.linalg.norm(ey)

    ez = np.cross(ex,ey)
    ez = ez/np.linalg.norm(ez)

    return np.asarray([ex,ey,ez]).T

def computeHelixMatrix(a,pitch,e,sigma):

    b     = pitch / (2.0 * np.pi)

    sOne = computeHelixParticlesDistance(a,b,e,sigma)

    R_0 = discreteHelixFrame(0,sOne,a,b,e)
    R_1 = discreteHelixFrame(1,sOne,a,b,e)

    return R_0.T@R_1 # R_1 in the basis of R_0

def energy(pos,ori,sigma,cNext,cPrevious):

    E = 0.0
    for i in range(N-1):
        pos_curr = pos[i]
        pos_next = pos[i+1]

        #connection_curr = ori[i]@np.asarray([sigma/2.0,0.0,0.0])+pos_curr
        #connection_next = ori[i+1]@np.asarray([-sigma/2.0,0.0,0.0])+pos_next # This is not correct, change in the future

        connection_curr = ori[i]@cNext+pos_curr
        connection_next = ori[i+1]@cPrevious+pos_next # This is not correct, change in the future

        r = np.linalg.norm(connection_curr-connection_next)

        Ebond = 0.5*Kb*r**2

        if np.abs(Ebond) < 1e-5:
            Ebond = 0.0

        E += Ebond

        ori_curr = ori[i]
        ori_next = ori[i+1]

        # We compute the expected orientation of the next particle given the current orientation
        ori_next_expected = ori_curr@R_H

        ex_next_expected = ori_next_expected[:,0]
        ey_next_expected = ori_next_expected[:,1]
        ez_next_expected = ori_next_expected[:,2]

        ex_next = ori_next[:,0]
        ey_next = ori_next[:,1]
        ez_next = ori_next[:,2]

        Ecurr_next = Ka*(1.0-np.dot(ex_next,ex_next_expected))+Kd*(1.0-np.dot(ey_next,ey_next_expected))

        # We compute the expected orientation of the current particle given the next orientation
        ori_curr_expected = ori_next@R_H.T

        ex_curr_expected = ori_curr_expected[:,0]
        ey_curr_expected = ori_curr_expected[:,1]
        ez_curr_expected = ori_curr_expected[:,2]

        ex_curr = ori_curr[:,0]
        ey_curr = ori_curr[:,1]
        ez_curr = ori_curr[:,2]

        Enext_curr = Ka*(1.0-np.dot(ex_curr,ex_curr_expected))+Kd*(1.0-np.dot(ey_curr,ey_curr_expected))

        # If Ecurr_next or Enext_curr are small enough, return 0.0
        if Ecurr_next < 1e-6:
            Ecurr_next = 0.0
        if Enext_curr < 1e-6:
            Enext_curr = 0.0

        E += Ecurr_next + Enext_curr

        #print("Ebonds:",Ebond,"Ecurr_next:",Ecurr_next,"Enext_curr:",Enext_curr)

    return E

def computeConnections(a,pitch,helicity,sigma):

    R_H = computeHelixMatrix(a,pitch,helicity,sigma)

    b     = pitch / (2.0 * np.pi)

    sOne = computeHelixParticlesDistance(a,b,helicity,sigma)

    R_0 = discreteHelixFrame(0,sOne,a,b,helicity)
    R_1 = discreteHelixFrame(1,sOne,a,b,helicity)

    pos0 = helixEquation(0.0,a,b,helicity)
    pos1 = helixEquation(sOne,a,b,helicity)

    pos0 = R_0.T@pos0
    pos1 = R_0.T@pos1

    dr   = (pos1-pos0)
    dr   = dr/np.linalg.norm(dr)
    dr   = dr*sigma/2.0

    connectionNext = dr

    pos0 = helixEquation(0.0,a,b,helicity)
    pos1 = helixEquation(sOne,a,b,helicity)

    pos0 = R_1.T@pos0
    pos1 = R_1.T@pos1

    dr   = (pos0-pos1)
    dr   = dr/np.linalg.norm(dr)
    dr   = dr*sigma/2.0

    connectionPrevious = dr

    return connectionNext,connectionPrevious

def generateHelix(N,a,pitch,helicity,sigma,initPos = np.array([0.0,0.0,0.0]),initOri = np.eye(3)):
    """
    Generates a helix with N beads, pitch and helicity
    """

    R_H = computeHelixMatrix(a,pitch,helicity,sigma)

    pos = np.zeros((N,3))
    ori = np.zeros((N,3,3))

    pos[0]     = initPos
    ori[0,:,:] = initOri

    for i in range(1,N):
        pos[i]  = ori[i-1]@np.asarray([sigma,0.0,0.0])+pos[i-1]
        ori[i]  = ori[i-1]@R_H

    return pos,ori

def helixAngularPerturbation(ori,angleRange):

    N = ori.shape[0]

    oriPerturbed = np.zeros((N,3,3))
    for i in range(len(oriPerturbed)):
        perturbation = smallRandomRotation(angleRange).as_matrix()
        oriPerturbed[i] = ori[i]@perturbation

    return oriPerturbed

def writeHelix(pos,ori,sigma,filename):

    with open(filename, "w") as f:
        f.write("#\n")
        for i in range(len(pos)):

            p  = pos[i,:]

            ex = ori[i][:,0]*sigma/2 + p
            ey = ori[i][:,1]*sigma/2 + p
            ez = ori[i][:,2]*sigma/2 + p

            f.write(f"{p[0] } {p[1] } {p[2] } {sigma/2.0 } 0 \n")
            f.write(f"{ex[0]} {ex[1]} {ex[2]} {sigma/10.0} 1 \n")
            f.write(f"{ey[0]} {ey[1]} {ey[2]} {sigma/10.0} 2 \n")
            f.write(f"{ez[0]} {ez[1]} {ez[2]} {sigma/10.0} 3 \n")

def writeHelixStream(pos,ori,sigma,fhandle):

    fhandle.write("#\n")

    for i in range(len(pos)):
        p  = pos[i,:]

        ex = ori[i][:,0]*sigma/2 + p
        ey = ori[i][:,1]*sigma/2 + p
        ez = ori[i][:,2]*sigma/2 + p

        fhandle.write(f"{p[0] } {p[1] } {p[2] } {sigma/2.0 } 0 \n")
        fhandle.write(f"{ex[0]} {ex[1]} {ex[2]} {sigma/10.0} 1 \n")
        fhandle.write(f"{ey[0]} {ey[1]} {ey[2]} {sigma/10.0} 2 \n")
        fhandle.write(f"{ez[0]} {ez[1]} {ez[2]} {sigma/10.0} 3 \n")
