Overview
========

UAMMD-structured is designed with the primary goal of simulating **generic physical systems**. This process is accomplished through the use of several entities, each responsible for managing and storing some of the information needed to run a simulation. Each of these entities is responsible for different aspects of the simulation, such as managing particle states, handling interactions and performing the integration process.

A key aspect of UAMMD-structured is its strong reliance on input. As highlighted earlier, there are two primary ways to operate UAMMD-structured. It can be compiled into a single, general-purpose binary that processes inputs in a specific format. Alternatively, it can be executed through a Python interface, where input is provided in the form of a dictionary. In both cases the format of the input is similar. The customization of the simulation is achieved by modifying the options available in the input file or dictionary. Through this input, the properties of the following basic entities are specified, determining the type of simulation being conducted:

- **System:** Specifies technical aspects of the simulation, like the name of the simulation, backup and restarting systems, or the seed for random number generation.

- **Integrator:** This section defines the integrators employed in the simulation. Integrators are algorithms responsible for updating the state of the simulation. UAMMD-structured provides a range of integrators, from the classical Verlet algorithm, to more complex ones used for simulating systems with fluid interactions. Commonly, only one integrator is used; however, UAMMD-structured also offers the option to sequentially apply multiple integrators.

- **Global:** Adds elements to the simulation not related to individual particles, such as the units used, the type of particles, or the ensemble (temperature, volume ...) for the simulation.

- **State:** Specifies the properties of particles that may vary throughout the simulation or that are unique to each particle, like position or orientation.

- **Topology:** Encapsulates information typically found in topology files in other software (e.g., GROMACS). This facilitates entry for new users by using a common language in the field. Topology specifies:

  - **Structure:** Associates each particle with a type and specific structures, like the residue or the chain it belongs to.

  - **Force Field:** Defines the interactions to which the particles are subjected, such as bonds, short-range interactions, external fields ...

- **Simulation Steps:** Selects the type of information obtained from the simulation. These are operations performed at user-defined intervals, involving different measurements or writing information to a file, such as dumping the state of the simulation or measuring temperature and pressure.

.. figure:: /img/simulation.png
   :alt: UAMMD-structured architecture

   This figure illustrates the architecture of the UAMMD-structured simulation framework. It highlights the key components involved in simulating physical systems. The architecture is broken down into several main entities: *System*, *Integrator*, *Global*, *State*, *Topology*, and *Simulation Steps*.

Upon specifying these entities, UAMMD-structured processes, connects them, and initiates the simulation. **Therefore, understanding UAMMD-structured is equivalent to knowing the different options available for these entities**.
