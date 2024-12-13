import sympy as sp
import numpy as np
import quaternion

import json

# Import lambdify
from sympy.utilities.lambdify import lambdify

# Define rotation matrices A,B and R
A = sp.Matrix([
    [sp.Symbol('a_1x'), sp.Symbol('a_2x'), sp.Symbol('a_3x')],
    [sp.Symbol('a_1y'), sp.Symbol('a_2y'), sp.Symbol('a_3y')],
    [sp.Symbol('a_1z'), sp.Symbol('a_2z'), sp.Symbol('a_3z')]
])

B = sp.Matrix([
    [sp.Symbol('b_1x'), sp.Symbol('b_2x'), sp.Symbol('b_3x')],
    [sp.Symbol('b_1y'), sp.Symbol('b_2y'), sp.Symbol('b_3y')],
    [sp.Symbol('b_1z'), sp.Symbol('b_2z'), sp.Symbol('b_3z')]
])

R = sp.Matrix([
    [sp.Symbol('R_1x'), sp.Symbol('R_2x'), sp.Symbol('R_3x')],
    [sp.Symbol('R_1y'), sp.Symbol('R_2y'), sp.Symbol('R_3y')],
    [sp.Symbol('R_1z'), sp.Symbol('R_2z'), sp.Symbol('R_3z')]
])

# Define K
K = sp.symbols('K')

def vee_operator(matrix):
    # Check if matrix is 3x3
    if matrix.shape != (3, 3):
        raise ValueError("Matrix must be 3x3")

    # Check if matrix is skew-symmetric: M = -M^T
    # For SymPy matrices, we need to check each element
    sum_matrix = (matrix + matrix.transpose())
    is_skew = all(sum_matrix[i,j].simplify() == 0
                  for i in range(3) for j in range(3))

    if not is_skew:
        raise ValueError("Matrix is not skew-symmetric")

    # Apply vee operator: extract (3,2), (1,3), (2,1) elements
    return sp.Matrix([
        matrix[2, 1],  # (3,2) element
        matrix[0, 2],  # (1,3) element
        matrix[1, 0]   # (2,1) element
    ])

def pi_func(matrix_A,matrix_B):
    return vee_operator(matrix_A@matrix_B.transpose()-matrix_B@matrix_A.transpose())

def frob(matrix_A,matrix_B):
    return (matrix_A.transpose()@matrix_B).trace()

def pot(A,B,R,K):
    RR = A.transpose()@B
    fp = frob(RR,R)
    return K*(3.0-fp)/4.0

def matrix_derivative(f, X):
    # Create a matrix of the same size as X to store derivatives
    rows, cols = X.shape
    dfdX = sp.zeros(rows, cols)

    # For each element in X, compute the partial derivative
    for i in range(rows):
        for j in range(cols):
            # Take derivative with respect to each element of X
            dfdX[i,j] = f.diff(X[i,j])

    return dfdX

energy = pot(A,B,R,K)
ta     = pi_func(A,matrix_derivative(energy,A))
tb     = pi_func(B,matrix_derivative(energy,B))

# Lambdify the expressions
var_list = [A[0,0], A[0,1], A[0,2],
            A[1,0], A[1,1], A[1,2],
            A[2,0], A[2,1], A[2,2],
            B[0,0], B[0,1], B[0,2],
            B[1,0], B[1,1], B[1,2],
            B[2,0], B[2,1], B[2,2],
            R[0,0], R[0,1], R[0,2],
            R[1,0], R[1,1], R[1,2],
            R[2,0], R[2,1], R[2,2],
            K]

print("Lambdifying expressions...")
energy_func = lambdify(var_list, energy, modules='numpy')
ta_func = lambdify(var_list, ta, modules='numpy')
tb_func = lambdify(var_list, tb, modules='numpy')

###############################################

paramFile = "../parameters.json"
stateFile = "output.itpd"

with open(paramFile) as f:
    params = json.load(f)

K_val = params['Krap']
R_q   = params['R']
R_q   = np.quaternion(R_q[0], R_q[1], R_q[2], R_q[3])

particles = []
skiplines = 2
with open(stateFile) as f:
    # Format is id type x y z q0 q1 q2 q3
    for line in f:
        if skiplines > 0:
            skiplines -= 1
            continue
        parts = line.split()

        id_ = int(parts[0])
        type_ = int(parts[1])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        q0 = float(parts[5])
        q1 = float(parts[6])
        q2 = float(parts[7])
        q3 = float(parts[8])

        particles.append(np.quaternion(q0, q1, q2, q3))

#Transform quaternion to rotation matrix
A_val = quaternion.as_rotation_matrix(particles[0])
B_val = quaternion.as_rotation_matrix(particles[1])
R_val = quaternion.as_rotation_matrix(R_q)

#############################################

func_args = [A_val[0,0], A_val[0,1], A_val[0,2],
             A_val[1,0], A_val[1,1], A_val[1,2],
             A_val[2,0], A_val[2,1], A_val[2,2],
             B_val[0,0], B_val[0,1], B_val[0,2],
             B_val[1,0], B_val[1,1], B_val[1,2],
             B_val[2,0], B_val[2,1], B_val[2,2],
             R_val[0,0], R_val[0,1], R_val[0,2],
             R_val[1,0], R_val[1,1], R_val[1,2],
             R_val[2,0], R_val[2,1], R_val[2,2],
             K_val]

energy_val = energy_func(*func_args)/2.0 #Divide by 2 because the energy is distributed between the two particles
ta_val = ta_func(*func_args)
tb_val = tb_func(*func_args)

print("Energy: ", energy_val)
print("Torque A: ", ta_val[0], ta_val[1], ta_val[2])
print("Torque B: ", tb_val[0], tb_val[1], tb_val[2])








