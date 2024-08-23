Simulation
==========

In UAMMD-structured, a simulation is defined through a structured input that specifies various components of the system and the simulation process. While we've seen a concrete example in the "First Example" section, let's break down the general structure and components of a UAMMD-structured simulation.

Key Components
--------------

A typical UAMMD-structured simulation input consists of the following main sections:

1. **System**:

   Specifies basic simulation parameters such as the name of the simulation.

2. **Global**:

   Sets up fundamental aspects of the simulation environment, including:
   - Unit system
   - Particle types and their properties
   - Ensemble parameters (e.g., box size, temperature)

3. **Integrator**:

   Defines the algorithm(s) used to evolve the system over time. This section includes:
   - Type of integrator (e.g., Langevin, Brownian dynamics)
   - Integration parameters (e.g., time step, friction constant)
   - Schedule for applying different integrators if more than one is used

4. **State**:

   Specifies the initial configuration of the system, including:
   - Particle positions
   - Velocities (optional)
   - Other particle-specific properties

5. **Topology**:

   Defines the structure of the system and the interactions between particles. It's divided into two main subsections:
   - **Structure**: Associates particles with types and higher-level structures (e.g., residues, chains)

   - **Force Field**: Specifies the interactions present in the system, such as:

     - Bonded interactions (e.g., harmonic bonds, angle potentials)
     - Non-bonded interactions (e.g., Lennard-Jones potential)
     - External fields or constraints

6. **Simulation Steps**:

   Defines periodic operations to be performed during the simulation, such as:
   - Outputting system information
   - Saving system state
   - Calculating and recording specific measurements

Structure and Formatting
------------------------

The simulation input is structured hierarchically, with each main section containing various subsections and entries. As explained in the Input System section, this structure is typically represented in YAML or JSON format, with each component defined by its type, parameters, and (where applicable) data.

It's important to note that the order of sections in the input file is not significant, thanks to the YAML/JSON format. Each entity resides in its unique section, and within these sections, there might be one or several data entries.

Every data entry (excluding State) possesses a "type" field, which is compulsory and must always be present. However, in some entries, we have both parameters and data, while in others, we might only have parameters or just data.

Flexibility and Customization
-----------------------------

This structure allows for great flexibility in defining simulations. Users can easily:

- Modify system parameters
- Change interaction potentials
- Adjust integration schemes
- Add or remove measurement and output operations

By adjusting these components, users can tailor the simulation to a wide variety of physical systems and research questions, from simple particle systems to complex biomolecular simulations.

Conclusion
----------

Understanding this structure is key to effectively using UAMMD-structured. By mastering the organization and options available in each section, users can construct sophisticated simulations tailored to their specific research needs. The example provided in the "First Example" section serves as a practical illustration of how these components come together to define a complete simulation.
