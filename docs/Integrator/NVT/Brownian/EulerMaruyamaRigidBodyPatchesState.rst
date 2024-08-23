EulerMaruyamaRigidBodyPatchesState
----------------------------------

The EulerMaruyamaRigidBodyPatchesState integrator extends the EulerMaruyamaRigidBody integrator to include support for patchy particles with state transitions. This integrator is designed for simulations of particles with directional interactions and internal states.

----

* **type**: ``Brownian``, ``EulerMaruyamaRigidBodyPatchesState``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the implicit solvent :math:`[mass/(distance \cdot time)]` (optional if translational and rotational self-diffusion are defined in particle data)

Example:

.. code-block::

   "eulerMaruyamaRigidBodyPatchesState":{
     "type":["Brownian","EulerMaruyamaRigidBodyPatchesState"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0
     }
   }

.. note::
   This integrator is specifically designed to work with patchy particle interactors that support state transitions. It updates both the positions and orientations of particles, as well as their internal states.

.. warning::
   This integrator requires that patchy particle interactors are properly set up and added to the simulation. It assumes that particle masses, radii, orientations (as quaternions), and patch states are defined in the particle data.

Functionality:
   1. Performs Brownian dynamics integration for particle positions and orientations.
   2. Updates the transition probabilities for patch states.
   3. Applies state transi
