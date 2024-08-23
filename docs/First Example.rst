First Example
=============

To demonstrate how UAMMD-structured works in practice, let's walk through a simple example simulation. This example will showcase a system of four particles configured in a simple rod-like structure, where consecutive particles are bonded by harmonic springs. Additionally, a harmonic angular potential is added between the first three and the last three particles.

Here's the complete input for this simulation:

.. code-block:: yaml

   system:
     info:
       type: ["Simulation", "Information"]
       parameters:
         name: "ROD"

   global:
     units:
       type: ["Units", "None"]
     types:
       type: ["Types", "Basic"]
       labels: ["name", "mass", "radius", "charge"]
       data:
         - ["A", 1.0, 0.5, 0.0]
     ensemble:
       type: ["Ensemble", "NVT"]
       labels: ["box", "temperature"]
       data:
         - [[[10.0, 10.0, 10.0], 1.0]]

   integrator:
     lang_bbk:
       type: ["Langevin", "BBK"]
       parameters:
         timeStep: 0.001
         frictionConstant: 1.0
     schedule:
       type: ["Schedule", "Integrator"]
       labels: ["order", "integrator", "steps"]
       data:
         - [1, "lang_bbk", 1000000]

   state:
     labels: ["id", "position"]
     data:
       - [0, [0.0, 0.0, 0.0]]
       - [1, [0.0, 0.0, 1.0]]
       - [2, [0.0, 0.0, 2.0]]
       - [3, [0.0, 0.0, 3.0]]

   topology:
     structure:
       labels: ["id", "type", "modelId"]
       data:
         - [0, "A", 0]
         - [1, "A", 0]
         - [2, "A", 0]
         - [3, "A", 0]
     forceField:
       bonds:
         type: ["Bond2", "HarmonicCommon_K_r0"]
         parameters:
           K: 100.0
           r0: 1.0
         labels: ["id_i", "id_j"]
         data:
           - [0, 1]
           - [1, 2]
           - [2, 3]
       angles:
         type: ["Bond3", "KratkyPorodCommon_K"]
         parameters:
           K: 50.0
         labels: ["id_i", "id_j", "id_k"]
         data:
           - [0, 1, 2]
           - [1, 2, 3]

   simulationStep:
     info:
       type: ["UtilsStep", "InfoStep"]
       parameters:
         intervalStep: 10000
     saveState:
       type: ["WriteStep", "WriteStep"]
       parameters:
         outputFilePath: "test"
         outputFormat: "sp"
         intervalStep: 10000

Let's break down this example to understand each component:

1. **System**: Specifies the name of the simulation as "ROD".

2. **Global**: Sets up the basic parameters of the simulation:
   
   - Uses no specific unit system
   - Defines a single particle type "A" with mass 1.0, radius 0.5, and no charge
   - Sets up an NVT ensemble with a cubic box of side length 10.0 and temperature 1.0

3. **Integrator**: Uses a Langevin BBK integrator with a time step of 0.001 and friction constant of 1.0. The simulation is scheduled to run for 1,000,000 steps.

4. **State**: Defines the initial positions of the four particles along the z-axis.

5. **Topology**: 

   - Structure: Assigns all particles to type "A" and model 0.
   - Force Field: 
     - Defines harmonic bonds between consecutive particles with spring constant 100.0 and equilibrium distance 1.0.
     - Adds Kratky-Porod angle potentials with constant 50.0 for the two sets of three consecutive particles.

6. **Simulation Step**: Sets up two periodic operations:

   - Prints simulation information every 10,000 steps
   - Saves the simulation state every 10,000 steps

This example demonstrates how UAMMD-structured uses a structured input to define all aspects of a molecular dynamics simulation, from system setup to runtime behavior.
