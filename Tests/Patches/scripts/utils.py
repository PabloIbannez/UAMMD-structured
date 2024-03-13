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
