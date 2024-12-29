import os

import json

import numpy as np
import quaternion

import matplotlib.pyplot as plt

from tqdm import tqdm

import subprocess

paramsFile   = "../parameters.json"
with open(paramsFile) as f:
    params = json.load(f)

temperature = params["temperature"]
beta = 1.0 / temperature

MC_steps = 10000000
MC_equilibration = 1000
MC_output = 100000

t_skip = 200

def frob_product(q1, q2, qR):
    q1_read = quaternion.from_float_array([q1[0], q1[1], q1[2], q1[3]])
    q2_read = quaternion.from_float_array([q2[0], q2[1], q2[2], q2[3]])

    qR_read = quaternion.from_float_array([qR[0], qR[1], qR[2], qR[3]])

    # Check if norm is 1
    assert np.isclose(q1_read.norm(), 1)
    assert np.isclose(q2_read.norm(), 1)
    assert np.isclose(qR_read.norm(), 1)

    A = quaternion.as_rotation_matrix(q1_read)
    B = quaternion.as_rotation_matrix(q2_read)

    R = quaternion.as_rotation_matrix(qR_read)

    RR = A.T @ B

    frob = (RR.T @ R).trace()

    return frob

def energy(q1, q2, qR, K):
    frob = frob_product(q1, q2, qR)
    E = K * (3 - frob) / 4
    return E

def frob_angle(q1, q2, qR):
    frob = frob_product(q1, q2, qR)
    theta = np.arccos((frob - 1) / 2)
    return theta


N  = params["N"]
qR = params["R"]
K  = params["Krap"]

# Generate outputMC.itpd.
# We execute MC program, the parameters are: N T K q0 q1 q2 q3 steps equilibration output filename
subprocess.run(["../scripts/MC", str(N),
                                 str(temperature), str(K),
                                 str(qR[0]), str(qR[1]), str(qR[2]), str(qR[3]),
                                 str(MC_steps), str(MC_equilibration), str(MC_output),
                                 "outputMC.itpd"])

files = ["output.itpd", "outputMC.itpd"]

for filepath in files:

    print(f"Processing {filepath}")

    # Remove all lines that start with N or are empty
    with open(filepath, "r") as f:
        with open("tmp.dat", "w") as f2:
            for line in f:
                if not line.startswith(str(N)) and not line.isspace():
                    f2.write(line)
    data = np.loadtxt("tmp.dat")
    # Remove temporary file
    os.remove("tmp.dat")

    # Data are lines of the form:
    # id type x y z qw qx qy qz (9 fields)
    # split them into chunks of size N
    data = data.reshape(-1, N, 9)

    # Remove the first t_skip frames
    data = data[t_skip:]

    theta_data = []
    for t,frame in enumerate(tqdm(data)):
        for i in range(N-1):
            q1 = frame[i, 5:]
            q2 = frame[i+1, 5:]

            theta = frob_angle(q1, q2, qR)
            theta_data.append(theta)

    # Plot the histogram of sampled angles
    plt.hist(theta_data, bins=100, density=True, alpha=0.5, label=f"{filepath}")

plt.xlim(0, np.pi)
plt.xlabel("Theta")
plt.ylabel("Probability Density")
plt.legend()
plt.show()
