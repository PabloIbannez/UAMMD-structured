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
# # Worm-Like Chain (WLC) Compression Simulation Tutorial for UAMMD-structured

# %% [markdown]
# ## Introduction
#
# This notebook demonstrates how to create and run a simulation of polymer chains under compression using the Worm-Like Chain (WLC) model in UAMMD-structured. We'll use a spherical shell that compresses the polymers over time. We'll walk through the process step-by-step, explaining each part of the simulation setup in detail.

# %% [markdown]
# ## Setup
#
# First, let's import the necessary libraries:

# %%
import os
import numpy as np
import pyUAMMD

# os: Used for directory operations
# numpy: Used for numerical operations and array manipulations
# pyUAMMD: The Python interface for UAMMD-structured

# %% [markdown]
# ## Simulation Parameters
#
# Here we define the key parameters for our WLC compression simulation:

# %%
# Number of beads per polymer chain
nBeads = 30

# Number of polymer chains
nPolymers = 100

# Bead diameter (and bond length)
sigma = 1.0

# Calculate box size (add some extra space)
L = nBeads*sigma + 10.0*sigma

# Simulation parameters
timeStep = 0.001
frictionConstant = 1.0

# Total number of simulation steps
nSteps = 100000*20000

# Frequency of information output
nStepsInfo = 1000

# Frequency of trajectory output
nStepsOutput = 10000

# Force constants for bonds and angles
Kb = 100.0  # Bond strength
Ka = 100.0  # Angle strength (bending rigidity)

# Compression velocity of the spherical shell
compressionVelocity = -0.01

print("Creating a WLC compression simulation with:")
print(f" - Number of beads per polymer: {nBeads}")
print(f" - Number of polymers: {nPolymers}")
print(f" - Bead diameter: {sigma}")
print(f" - Box size: {L}")
print(f" - Time step: {timeStep}")
print(f" - Friction constant: {frictionConstant}")
print(f" - Total number of steps: {nSteps}")
print(f" - Bond strength (Kb): {Kb}")
print(f" - Bending rigidity (Ka): {Ka}")
print(f" - Compression velocity: {compressionVelocity}")
print(f" - Output information every {nStepsInfo} steps")
print(f" - Output trajectory every {nStepsOutput} steps")

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
        "parameters": {"name": "WormLikeChainCompression"}
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
        "data": [["A", 1.0, sigma/2.0, 0.0]]
    },

    # Set the ensemble (NVT: constant Number of particles, Volume, and Temperature)
    "ensemble": {
        "type": ["Ensemble", "NVT"],
        "labels": ["box", "temperature"],
        "data": [[[L, L, L], 1.0]]
    }
}

# Set up the integrator (Langevin dynamics)
simulation["integrator"] = {
    "bbk": {
        "type": ["Langevin", "BBK"],
        "parameters": {
            "timeStep": timeStep,
            "frictionConstant": frictionConstant
        }
    },
    # Define the integration schedule
    "schedule": {
        "type": ["Schedule", "Integrator"],
        "labels": ["order", "integrator", "steps"],
        "data": [[1, "bbk", nSteps]]
    }
}

# %% [markdown]
# ## Initialize Particle Positions and Topology
#
# We'll place the polymer chains in the simulation box and set up the topology:

# %%
simulation["state"] = {
    "labels": ["id", "position"],
    "data": []
}

simulation["topology"] = {
    "structure": {
        "labels": ["id", "type", "modelId"],
        "data": []
    }
}

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads)):
        # Place beads along the z-axis, centered at the origin
        simulation["state"]["data"].append([particleId, [0.0, 0.0, sigma*j - nBeads*sigma/2.0]])
        simulation["topology"]["structure"]["data"].append([particleId, "A", i])
        particleId += 1

# %% [markdown]
# ## Define Force Field
#
# Set up the bonded interactions for our WLC model and the compressing spherical shell:

# %%
simulation["topology"]["forceField"] = {}

# Set up harmonic bonds
simulation["topology"]["forceField"]["bonds"] = {
    "type": ["Bond2", "Harmonic"],
    "parameters": {},
    "labels": ["id_i", "id_j", "K", "r0"],
    "data": []
}

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads) - 1):
        simulation["topology"]["forceField"]["bonds"]["data"].append([particleId, particleId + 1, Kb, sigma])
        particleId += 1
    particleId += 1

# Set up angle interactions (Kratky-Porod potential for bending rigidity)
simulation["topology"]["forceField"]["angles"] = {
    "type": ["Bond3", "KratkyPorod"],
    "parameters": {},
    "labels": ["id_i", "id_j", "id_k", "K"],
    "data": []
}

particleId = 0
for i in range(int(nPolymers)):
    for j in range(int(nBeads) - 2):
        simulation["topology"]["forceField"]["angles"]["data"].append([particleId, particleId + 1, particleId + 2, Ka])
        particleId += 1
    particleId += 2

# Set up the compressing spherical shell
simulation["topology"]["forceField"]["shell"] = {
    "type": ["External", "SphericalShell"],
    "parameters": {
        "shellCenter": [0.0, 0.0, 0.0],
        "shellRadius": nBeads*sigma + 2.0*sigma,
        "shellEpsilon": 1.0,
        "shellSigma": sigma,
        "minShellRadius": nBeads*sigma/2.0/2.0,
        "radiusVelocity": compressionVelocity
    }
}

# %% [markdown]
# ## Configure Simulation Steps
#
# Define what operations to perform during the simulation:

# %%
simulation["simulationStep"] = {
    # Output simulation information periodically
    "info": {
        "type": ["UtilsStep", "InfoStep"],
        "parameters": {"intervalStep": nStepsInfo}
    },
    # Save trajectory data periodically
    "output": {
        "type": ["WriteStep", "WriteStep"],
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
# This tutorial demonstrated how to set up and prepare a simulation of Worm-Like Chains under compression using UAMMD-structured. We covered:
# 1. Setting up the environment for Jupyter/Google Colab
# 2. Defining simulation parameters, including polymer-specific parameters and compression velocity
# 3. Creating the simulation object
# 4. Setting up the system, global parameters, and integrator
# 5. Initializing particle positions and topology for multiple polymer chains
# 6. Configuring the simulation topology and force field, including bonded and angle interactions
# 7. Setting up the compressing spherical shell
# 8. Setting up simulation steps for output
# 9. Writing the simulation file
# 10. Running the simulation
#
# Next steps could include:
# - Analyzing the output trajectory
# - Visualizing the polymer chains under compression
# - Calculating polymer properties such as end-to-end distance or radius of gyration as a function of compression
# - Modifying the simulation parameters to explore different polymer lengths, stiffnesses, or compression rates
# - Investigating the behavior of the system under different environmental conditions
