import numpy as np
import matplotlib.pyplot as plt
import pyUAMMD

plotData = False

data = np.loadtxt('low_resolution.dat')
Npx = 12
Npy = 13

# Check if data has size Npx*Npy
if data.shape[0] != Npx*Npy:
    raise ValueError("Data does not have the correct size")

if plotData:

    # Extract unique x and y values
    x_values = np.unique(data[:, 0])
    y_values = np.unique(data[:, 1])
    z_values = np.unique(data[:, 2])

    # Create a grid to reshape values
    X, Y = np.meshgrid(x_values, y_values)
    Z = np.zeros_like(X, dtype=float)

    # Fill Z with corresponding values
    for i in range(data.shape[0]):
        xi = np.where(x_values == data[i, 0])[0][0]
        yi = np.where(y_values == data[i, 1])[0][0]
        Z[yi, xi] = data[i, 2]

    # Create the heatmap
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
    plt.colorbar(label="Value")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Heatmap of Data")

    plt.show()

N  = 1000
Nb = 10

sigma = 0.5
gamma = 10.0
tipRadius = 25.0

L = 100.0

timeStep = 0.001
frictionConstant = 1.0

nSteps       = 100000
nStepsInfo   = 1000
nStepsOutput = 1000

simulation = pyUAMMD.simulation()

simulation["system"] = {
    "info": {
        "type": ["Simulation", "Information"],
        "parameters": {"name": "afmImage"}
    }
}

simulation["global"] = {
    "units": {"type": ["Units", "None"]},
    "types": {
        "type": ["Types", "Basic"],
        "labels": ["name", "mass", "radius", "charge"],
        "data": [["A", 1.0, 0.5, 0.0]]
    },
    "ensemble": {
        "type": ["Ensemble", "NVT"],
        "labels": ["box", "temperature"],
        "data": [[[L, L, L], 1.0]]
    }
}

simulation["integrator"] = {
    "bbk": {
        "type": ["Langevin", "BBK"],
        "parameters": {
            "timeStep": timeStep,
            "frictionConstant": frictionConstant
        }
    },
    "schedule": {
        "type": ["Schedule", "Integrator"],
        "labels": ["order", "integrator", "steps"],
        "data": [[1, "bbk", nSteps]]
    }
}

simulation["state"] = {
    "labels": ["id", "position"],
    "data": []
}

particleId = 0
for nb in range(Nb):
    for n in range(N):
        # Generate a random position in x,y in the rage [-L/2, L/2]
        x = np.random.uniform(-L/2.0, L/2.0)
        y = np.random.uniform(-L/2.0, L/2.0)
        simulation["state"]["data"].append([particleId, [x, y, 0.0]])
        particleId += 1

simulation["topology"] = {
    "structure": {
        "labels": ["id", "type", "batchId"],
        "data": []
    }
}

particleId = 0
for nb in range(Nb):
    for n in range(N):
        simulation["topology"]["structure"]["data"].append([particleId, "A",nb])
        particleId += 1

simulation["topology"]["forceField"] = {}

simulation["topology"]["forceField"]["afmimage"] = {
    "type": ["AFMimage", "Takada"],
    "parameters": {
        "Npx": Npx,
        "Npy": Npy,
        "sigma": sigma,
        "gamma": gamma,
        "tipRadius": tipRadius
    },
    "labels": ["x", "y", "height", "batchId"],
    "data": []
}

for nb in range(Nb):
    for pixel in range(Npx*Npy):
        x, y, z = data[pixel]
        simulation["topology"]["forceField"]["afmimage"]["data"].append([x, y, z, nb])


simulation["simulationStep"] = {
    "info": {
        "type": ["UtilsStep", "InfoStep"],
        "parameters": {"intervalStep": nStepsInfo}
    }
}

simulation.write("simulation.json")

