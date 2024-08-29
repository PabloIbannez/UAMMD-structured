# %%
# Seting up the environment if running in Jupyter
import sys
in_jupyter = 'ipykernel' in sys.modules

if in_jupyter:
    get_ipython().system('pip install pyUAMMD')

    get_ipython().system('apt install libomp-dev')

    get_ipython().system('gdown 1rQYRMRFAEdmv8UCN5gcIpJL95fCtIBJg')
    get_ipython().system('dpkg -i uammdstructured-1.0.0-Linux.deb')
else:
    # Not in Jupyter, so we assume that the library is already installed
    pass

# %% [markdown]
# # Ideal Gas Simulation Tutorial for UAMMD-structured

# %% [markdown]
# ## Introduction
#
# This notebook demonstrates how to create and run a simulation of an ideal gas using UAMMD-structured. We'll walk through the process step-by-step, explaining each part of the simulation setup in detail.

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
# Here we define the key parameters for our ideal gas simulation:

# %%
# Number of particles in the simulation
N = 2e5

# Concentration of particles in the box (particles per unit volume)
conc = 0.1

# Calculate box length based on number of particles and concentration
# We use the cube root to get the length of a cubic box
L = (N / conc)**(1.0/3.0)

# Time step for the simulation (in reduced units)
timeStep = 0.001

# Friction constant for the Langevin integrator
# This determines the strength of coupling to the heat bath
frictionConstant = 1.0

# Total number of simulation steps
# This determines the total simulation time: total_time = nSteps * timeStep
nSteps = 10000

# Frequency of information output (every nStepsInfo steps)
nStepsInfo = 1000

# Frequency of trajectory output (every nStepsOutput steps)
nStepsOutput = 10000

print("Creating a simulation with:")
print(f" - Number of particles: {N}")
print(f" - Concentration: {conc}")
print(f" - Box length: {L}")
print(f" - Time step: {timeStep}")
print(f" - Friction constant: {frictionConstant}")
print(f" - Total number of steps: {nSteps}")
print(f" - Output information every {nStepsInfo} steps")
print(f" - Output trajectory every {nStepsOutput} steps")
print()

# %% [markdown]
# ## Creating the Simulation
#
# Now, let's create our simulation object and set up its various components:

# %%
# Initialize the simulation object
# This object will contain all the information about our simulation
simulation = pyUAMMD.simulation()

# Set up the system information
# This includes basic metadata about the simulation
simulation["system"] = {
    "info": {
        "type": ["Simulation", "Information"],
        "parameters": {"name": "IdealGas"}
    }
}

# Define global parameters
# This includes units, particle types, and ensemble information
simulation["global"] = {
    # Set the unit system (in this case, we're using reduced units)
    "units": {"type": ["Units", "None"]},

    # Define particle types (here we have one type "A" with mass 1.0, radius 0.5, and no charge)
    "types": {
        "type": ["Types", "Basic"],
        "labels": ["name", "mass", "radius", "charge"],
        "data": [["A", 1.0, 0.5, 0.0]]
    },

    # Set the ensemble (NVT: constant Number of particles, Volume, and Temperature)
    # We define the box size and temperature (set to 1.0 in reduced units)
    "ensemble": {
        "type": ["Ensemble", "NVT"],
        "labels": ["box", "temperature"],
        "data": [[[L, L, L], 1.0]]
    }
}

# Set up the integrator
# We're using a Langevin integrator (BBK algorithm) for NVT ensemble
simulation["integrator"] = {
    "bbk": {
        "type": ["Langevin", "BBK"],
        "parameters": {
            "timeStep": timeStep,
            "frictionConstant": frictionConstant
        }
    },
    # Define the integration schedule (which integrator to use and for how many steps)
    "schedule": {
        "type": ["Schedule", "Integrator"],
        "labels": ["order", "integrator", "steps"],
        "data": [[1, "bbk", nSteps]]
    }
}

# Initialize particle positions
# We'll place particles on a grid within the simulation box
simulation["state"] = {
    "labels": ["id", "position"],
    "data": []
}

# Create a grid of positions within the box
Linterval = np.linspace(-(L-2.0)/2.0, (L-2.0)/2.0, int(np.ceil(np.cbrt(N))+1))

particleId = 0
for lx in Linterval:
    for ly in Linterval:
        for lz in Linterval:
            if particleId < N:
                simulation["state"]["data"].append([particleId, [lx, ly, lz]])
                particleId += 1

# Set up topology
# This defines the structure of our system (particle IDs and types)
simulation["topology"] = {
    "structure": {
        "labels": ["id", "type"],
        "data": [[i, "A"] for i in range(int(N))]
    }
}

# Configure simulation steps
# These define what operations to perform during the simulation
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
            "outputFormat": "dcd"
        }
    }
}

# %% [markdown]
# ## Writing the Simulation File
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
# This tutorial demonstrated how to set up and prepare a simulation of an ideal gas using UAMMD-structured. We covered:
# 1. Defining simulation parameters
# 2. Creating the simulation object
# 3. Setting up the system, global parameters, and integrator
# 4. Initializing particle positions
# 5. Configuring the simulation topology
# 6. Setting up simulation steps for output
# 7. Writing the simulation file
#
# Next steps could include:
# - Analyzing the output trajectory
# - Visualizing the results
# - Modifying the simulation parameters to explore different conditions
# - Adding interactions between particles to simulate a non-ideal gas

