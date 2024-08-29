Integrators
===========

Once the initial state of the simulation has been established through *State*, the next step is to determine how this state evolves. This is accomplished using integrators. In UAMMD-structured, integrators are the elements responsible for updating the state of the particles. Although other algorithms, like Monte Carlo, can be implemented, our focus here is on molecular dynamics integrators. Generally, these integrators evolve the state of the system over time by computing properties such as forces, which are then used to determine the position of the particles at subsequent simulation time points.

The selection of the integrator(s) and the number of integration steps required, are specified by the type of input given in the following example:

.. code-block:: yaml

   integrator:
     bbk1:
       parameters:
         frictionConstant: 1.0
         timeStep: 0.00001
       type: ["Langevin", "BBK"]
     bbk2:
       parameters:
         frictionConstant: 1.0
         timeStep: 0.001
       type: ["Langevin", "BBK"]
     schedule:
       data:
         - [1, "bbk1", 10000]
         - [2, "bbk2", 1000000]
       labels: ["order", "integrator", "steps"]
       type: ["Schedule", "Integrator"]

The integrators section must include an additional Data Entry named "**schedule**" of type "**Integrator**". This entry sets the sequence and the number of steps for each integrator. The integrators of the simulation are first defined (e.g., two Langevin integrators of the BBK subtype in the current example), followed by specifying their execution order and the number of steps in the schedule.

This ordered execution is particularly beneficial for certain processes such as equilibration. For example, in thermalization, an integrator with a smaller "dt" might be used initially for a few steps, followed by a longer period using an integrator with a larger "dt". This approach is also advantageous in fluid dynamics scenarios, where a Langevin integrator can first equilibrate the system, followed by a fluid integrator for analyzing dynamic properties.

UAMMD-structured implements integrators commonly used in coarse-grained model simulations, including a Langevin dynamics integrator, specifically one that implements the BBK method. It also includes Brownian dynamics integrators, including the Brownian dynamics integrator for rigid bodies. This integrator is essential for simulations of oriented particles such as in lipid models or patchy particles commonly used in coarse-grained modeling. Additionally, UAMMD-structured includes integrators not intended for molecular dynamics but used for generating initial positions. Sometimes when a system is generated, certain particles are placed in such a way that they cause a steric clash. To resolve this, we can employ certain algorithms that eliminate these clashes by slightly modifying the positions of the particles, particularly UAMMD-structured implements the steepest descent which places the system in a local energy minimum.

UAMMD stands out with its wide array of integrators, addressing various simulation levels from the microscopic to the hydrodynamic. At the microscopic level, the dynamics are governed by Newtonian equations of motion, with Molecular Dynamics (MD) simulations implemented using the velocity Verlet (VV) scheme. This level deals with the positions and momenta of all particles, covering scales of up to hundreds of nanometers and nanoseconds.

For Langevin dynamics, which represents the dynamics of solute particles in a thermal bath, UAMMD incorporates integrators like the Gronbech-Jensen integrator, Dissipative Particle Dynamics (DPD), and Smoothed Particle Hydrodynamics (SPH).

The hydrodynamic level in UAMMD replaces individual solvent particles with fluctuating hydrodynamic fields, using the Immersed Boundary (IB) method for fluctuating hydrodynamics, effective for minimal resolution particles (blobs). UAMMD offers both compressible and incompressible hydrodynamic schemes, the latter utilizing pseudo-spectral schemes and GPU-adapted Fast Fourier Transforms (FFT). At the Brownian level, UAMMD deals with the dynamics of solute particles with small inertia, solving Brownian Dynamics (BD) with and without hydrodynamic interactions. It includes various tools for handling complex operations associated with these simulations, like the Cholesky and Lanczos methods for calculating the noise-related square mobility matrix. Finally, UAMMD also incorporates integrators for specific boundary conditions, such as Quasi2D for doubly periodic slices of fluid and double periodic Stokes (DPStokes) for integrating hydrodynamic interactions in slabs. These integrations allow for efficient simulation of complex systems under varied boundary conditions.

All these integrators can be used through UAMMD-structured, which for these types of integrators acts merely as a wrapper to certain functions of UAMMD. They can be added in the integrators section, the types and parameters can be consulted in the online documentation. Furthermore, the simplicity of using these integrators, combined with the extensive range of interactions that can be added, makes UAMMD-structured an immensely useful tool for research in soft matter, enabling researchers to conduct complex and detailed simulations with ease.
