# %%
# Setting up the environment if running in Jupyter
import sys
in_jupyter = 'ipykernel' in sys.modules

if in_jupyter:
    # Install necessary packages
    get_ipython().system('pip install pyUAMMD')
    get_ipython().system('apt install libomp-dev')

    # Download and install UAMMD-structured
    get_ipython().system('gdown 1rQYRMRFAEdmv8UCN5gcIpJL95fCtIBJg')
    get_ipython().system('dpkg -i uammdstructured-1.0.0-Linux.deb')
else:
    # Not in Jupyter, so we assume that the library is already installed
    pass

# %% [markdown]
# # Patchy Particles Simulation Tutorial for UAMMD-structured

# %% [markdown]
# ## Introduction
#
# This notebook demonstrates how to create and run a simulation of patchy particles using UAMMD-structured. Patchy particles are particles with specific interaction sites (patches) on their surface, which can lead to interesting self-assembly behavior. We'll walk through the process step-by-step, explaining each part of the simulation setup in detail.

# %% [markdown]
# ## Setup
#
# First, let's import the necessary libraries:

# %%
import sys
import os
import json
from tqdm import tqdm
import numpy as np
import pyUAMMD
from scipy.spatial.transform import Rotation

# sys, os: Used for system and directory operations
# json: Used for JSON file operations
# tqdm: Used for progress bars
# numpy: Used for numerical operations and array manipulations
# pyUAMMD: The Python interface for UAMMD-structured
# scipy.spatial.transform: Used for generating random rotations

# %% [markdown]
# ## Simulation Parameters
#
# Here we define the key parameters for our patchy particles simulation:

# %%
# Number of particles
N = 2000

# Concentration of particles in the box
c = 0.3

# Particle diameter
sigma = 1.0

# Patch interaction parameters
E = 10.0  # Energy scale
rc = 0.5  # Cutoff radius
Kb = 1.1  # Spring constant

# Simulation parameters
dt = 0.001  # Time step
nSteps = 100000  # Total number of steps
nStepsOutput = 10000  # Frequency of trajectory output

# Avoid periodic boundary conditions
avoidPBC = False

# Define patch positions
connectionXup   = [ sigma/2.0,  0,  0]
connectionXdown = [-sigma/2.0,  0,  0]
connectionYup   = [ 0,  sigma/2.0,  0]
connectionYdown = [ 0, -sigma/2.0,  0]
connectionZup   = [ 0,  0,  sigma/2.0]
connectionZdown = [ 0,  0, -sigma/2.0]

print('Generating patchy particles simulation with parameters:')
print(f'N = {N}')
print(f'c = {c}')
print(f'sigma = {sigma}')
print(f'E = {E}')
print(f'rc = {rc}')
print(f'Kb = {Kb}')

# Calculate box size
L = (N/c)**(1.0/3.0)
box = [L, L, L]

if avoidPBC:
    box = [2.0*L, 2.0*L, 2.0*L]

# %% [markdown]
# ## Creating the Simulation
#
# Now, let's create our simulation object and set up its various components:

# %%
# Initialize the simulation object
simulation = pyUAMMD.simulation()

# Set up the system information
simulation["system"] = {
    "info": {
        "type": ["Simulation", "Information"],
        "parameters": {"name": "dynamicExponentialTest"}
    }
}

# Define global parameters
simulation["global"] = {
    # Set the unit system (in this case, we're using reduced units)
    "units": {"type": ["Units", "None"]},

    # Define particle types
    "types": {
        "type": ["Types", "Basic"],
        "labels": ["name", "mass", "radius", "charge"],
        "data": [["A", 1.0, sigma/2, 0.0]]
    },

    # Set the ensemble (NVT: constant Number of particles, Volume, and Temperature)
    "ensemble": {
        "type": ["Ensemble", "NVT"],
        "labels": ["box", "temperature"],
        "data": [[box, 1.0]]
    }
}

# Set up the integrator
simulation["integrator"] = {}

if dt == 0.0:
    simulation["integrator"]["integrator"] = {
        "type": ["None", "None"]
    }
else:
    simulation["integrator"]["integrator"] = {
        "type": ["Brownian", "EulerMaruyamaRigidBody"],
        "parameters": {
            "timeStep": dt,
            "viscosity": 1.0
        }
    }

simulation["integrator"]["schedule"] = {
    "type": ["Schedule", "Integrator"],
    "labels": ["order", "integrator", "steps"],
    "data": [[1, "integrator", nSteps]]
}

# %% [markdown]
# ## Initialize Particle Positions and Orientations
#
# We'll place particles randomly within the simulation box, ensuring no overlaps:

# %%
simulation["state"] = {
    "labels": ["id", "position", "direction"],
    "data": []
}

print("Inserting particles...")
for i in tqdm(range(N)):
    # Generate random position within the box
    clash = True
    while clash:
        x = np.random.uniform(-L/2, L/2)
        y = np.random.uniform(-L/2, L/2)
        z = np.random.uniform(-L/2, L/2)
        p = [x, y, z]

        if i == 0:
            clash = False

        for j in range(i):
            dx = p[0] - simulation["state"]["data"][j][1][0]
            dy = p[1] - simulation["state"]["data"][j][1][1]
            dz = p[2] - simulation["state"]["data"][j][1][2]

            # Consider PBC
            dx = dx - box[0]*np.round(dx/box[0])
            dy = dy - box[1]*np.round(dy/box[1])
            dz = dz - box[2]*np.round(dz/box[2])

            r2 = dx*dx + dy*dy + dz*dz

            if r2 < 0.8*sigma*sigma:
                clash = True
                break
            else:
                clash = False

    # Generate random orientation
    q = Rotation.random().as_quat()  # scalar last
    q = [q[3], q[0], q[1], q[2]]  # scalar first

    simulation["state"]["data"].append([i, p, q])

# %% [markdown]
# ## Set up Topology
#
# Define the structure of our system (particle IDs and types):

# %%
simulation["topology"] = {
    "structure": {
        "labels": ["id", "type"],
        "data": [[i, "A"] for i in range(N)]
    }
}

# %% [markdown]
# ## Define Force Field
#
# Set up the interactions for our patchy particles:

# %%
simulation["topology"]["forceField"] = {}

# Set up neighbor list
simulation["topology"]["forceField"]["nl"] = {
    "type": ["VerletConditionalListSet", "all"],
    "parameters": {"cutOffVerletFactor": 1.2}
}

# Set up WCA (Weeks-Chandler-Andersen) potential
simulation["topology"]["forceField"]["wca"] = {
    "type": ["NonBonded", "WCAType2"],
    "parameters": {
        "cutOffFactor": 2.5,
        "condition": "all"
    },
    "labels": ["name_i", "name_j", "epsilon", "sigma"],
    "data": [["A", "A", 1.0, sigma]]
}

# Set up patchy interactions
simulation["topology"]["forceField"]["exponential"] = {
    "type": ["PatchyParticles", "DynamicallyBondedPatchyParticles"],
    "patchesState": {
        "labels": ["id", "position"],
        "data": []
    },
    "patchesGlobal": {
        "fundamental": {
            "type": ["Fundamental", "DynamicallyBondedPatchyParticles"],
            "parameters": {"energyThreshold": 0.0}
        },
        "types": {
            "type": ["Types", "Basic"],
            "labels": ["name", "mass", "radius", "charge"],
            "data": [["P", 0.0, sigma/10.0, 0.0]]
        }
    },
    "patchesTopology": {
        "structure": {
            "labels": ["id", "type", "parentId"],
            "data": []
        },
        "forceField": {
            "verletList": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2}
            },
            "exponential": {
                "type": ["NonBondedPatches", "DistanceSwitchExponential"],
                "parameters": {"condition": "all"},
                "labels": ["name_i", "name_j", "E", "K", "rc"],
                "data": [["P", "P", E, Kb, rc]]
            }
        }
    }
}

# Add patches to each particle
for i in range(N):
    index = i * 6
    patches = [connectionXup, connectionXdown, connectionYup, connectionYdown, connectionZup, connectionZdown]
    for j, patch in enumerate(patches):
        patch_id = index + j
        simulation["topology"]["forceField"]["exponential"]["patchesState"]["data"].append([patch_id, patch])
        simulation["topology"]["forceField"]["exponential"]["patchesTopology"]["structure"]["data"].append([patch_id, "P", i])

# %% [markdown]
# ## Configure Simulation Steps
#
# Define what operations to perform during the simulation:

# %%
simulation["simulationStep"] = {
    # Output simulation information periodically
    "info": {
        "type": ["UtilsStep", "InfoStep"],
        "parameters": {"intervalStep": nStepsOutput}
    },
    # Save trajectory data periodically
    "output": {
        "type": ["WriteStep", "WritePatchyParticlesStep"],
        "parameters": {
            "intervalStep": nStepsOutput,
            "outputFilePath": "output",
            "outputFormat": "sp"
        }
    }
}

# %% [markdown]
# ## Write the Simulation File
#
# Finally, let's write our simulation to a JSON file:

# %%
print()
print("Writing simulation file...")
simulation.write("simulation.json")
print("Simulation file created successfully!")
print()

# %% [markdown]
# ## Running the Simulation
#
# To run the simulation, you would typically use the UAMMD-structured executable with the generated JSON file.

# %%
if in_jupyter:
    get_ipython().system('UAMMDlauncher simulation.json')
else:
    print()
    print("You can run the code now using: UAMMDlauncher simulation.json")

# %% [markdown]
# ## Conclusion
#
# This tutorial demonstrated how to set up and prepare a simulation of patchy particles using UAMMD-structured. We covered:
# 1. Setting up the environment for Jupyter/Google Colab
# 2. Defining simulation parameters, including patch-specific parameters
# 3. Creating the simulation object
# 4. Setting up the system, global parameters, and integrator
# 5. Initializing particle positions and orientations
# 6. Configuring the simulation topology and force field, including patchy interactions
# 7. Setting up simulation steps for output
# 8. Writing the simulation file
# 9. Running the simulation
#
# Next steps could include:
# - Analyzing the output trajectory
# - Visualizing the results, particularly the formation of structures due to patchy interactions
# - Modifying the simulation parameters to explore different patch configurations or interaction strengths
# - Investigating the self-assembly behavior of the patchy particles under various conditions
