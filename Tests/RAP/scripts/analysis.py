import numpy as np
import quaternion

import json

paramsFile   = "../parameters.json"
with open(paramsFile) as f:
    params = json.load(f)

def frob(q1, q2, R):
    q1_read = quaternion.from_float_array([q1[3], q1[0], q1[1], q1[2]])
    q2_read = quaternion.from_float_array([q2[3], q2[0], q2[1], q2[2]])

    qR = quaternion.from_float_array([R[3], R[0], R[1], R[2]])

    # Check if norm is 1
    assert np.isclose(q1_read.norm(), 1)
    assert np.isclose(q2_read.norm(), 1)
    assert np.isclose(qR.norm(), 1)

    A = quaternion.as_rotation_matrix(q1_read)
    B = quaternion.as_rotation_matrix(q2_read)

    RR = quaternion.as_rotation_matrix(qR)

    print(A,B,RR)

    return 0


N = params["N"]
R = params["R"]

filepath = "output.itpd"

# Remove all lines that start with N or are empty
with open(filepath, "r") as f:
    with open("tmp.dat", "w") as f2:
        for line in f:
            if not line.startswith(str(N)) and not line.isspace():
                f2.write(line)

data = np.loadtxt("tmp.dat")

# Data are lines of the form:
# id type x y z qx qy qz qw (9 fields)
# split them into chunks of size N
data = data.reshape(-1, N, 9)

for frame in data:
    for i in range(N-1):
        q1 = frame[i, 5:]
        q2 = frame[i+1, 5:]

        frob_prod = frob(q1, q2, R)

    input()







