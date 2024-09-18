# %% [markdown]
# # Lennard-Jones Fluid Simulation Tutorial for UAMMD-structured

# %% [markdown]
# ## Introduction
#
# This notebook demonstrates how to create and run a simulation of a Lennard-Jones fluid using UAMMD-structured. We'll walk through the process step-by-step, explaining each part of the simulation setup in detail.

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
# Here we define the key parameters for our Lennard-Jones fluid simulation:

# %%
# Number of particles in the simulation
N = 2e5

# Concentration of particles in the box (particles per unit volume)
conc = 0.1

# Calculate box length based on number of particles and concentration
L = (N / conc)**(1.0/3.0)

# Lennard-Jones particle diameter
sigma = 1.0

# Time step for the simulation
timeStep = 0.001

# Friction constant for the Langevin integrator
frictionConstant = 1.0

# Total number of simulation steps
nSteps = 100000*20000

# Frequency of information output (every nStepsInfo steps)
nStepsInfo = 1000

# Frequency of trajectory output (every nStepsOutput steps)
nStepsOutput = 10000

print("Creating a simulation with:")
print(f" - Number of particles: {N}")
print(f" - Concentration: {conc}")
print(f" - Box length: {L}")
print(f" - Particle diameter: {sigma}")
print(f" - Time step: {timeStep}")
print(f" - Friction constant: {frictionConstant}")
print(f" - Total number of steps: {nSteps}")
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
        "parameters": {"name": "LennardJones"}
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
# ## Initialize Particle Positions
#
# We'll place particles on a grid within the simulation box:

# %%
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

# %% [markdown]
# ## Set up Topology
#
# Define the structure of our system (particle IDs and types):

# %%
simulation["topology"] = {
    "structure": {
        "labels": ["id", "type"],
        "data": [[i, "A"] for i in range(int(N))]
    }
}

# %% [markdown]
# ## Define Force Field
#
# Set up the Lennard-Jones interaction and neighbor list:

# %%
simulation["topology"]["forceField"] = {}

# Set up neighbor list
simulation["topology"]["forceField"]["nl"] = {
    "type": ["VerletConditionalListSet", "all"],
    "parameters": {"cutOffVerletFactor": 1.2}
}

# Set up Lennard-Jones interaction
simulation["topology"]["forceField"]["lj"] = {
    "type": ["NonBonded", "LennardJonesType2"],
    "parameters": {
        "condition": "all",
        "cutOffFactor": 2.5
    },
    "labels": ["name_i", "name_j", "epsilon", "sigma"],
    "data": [["A", "A", 1.0, sigma]]
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
print()
print("You can run the code now using: UAMMDlauncher simulation.json")

# %% [markdown]
# ## Conclusion
#
# This tutorial demonstrated how to set up and prepare a simulation of a Lennard-Jones fluid using UAMMD-structured. We covered:
# 1. Defining simulation parameters
# 3. Creating the simulation object
# 4. Setting up the system, global parameters, and integrator
# 5. Initializing particle positions
# 6. Configuring the simulation topology and force field
# 7. Setting up simulation steps for output
# 8. Writing the simulation file
# 9. Running the simulation
#
# Next steps could include:
# - Analyzing the output trajectory
# - Visualizing the results
# - Modifying the simulation parameters to explore different conditions
