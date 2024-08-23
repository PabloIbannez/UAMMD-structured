EulerMaruyamaRigidBody
----------------------

The EulerMaruyamaRigidBody integrator extends the Euler-Maruyama method to simulate rigid body Brownian dynamics [Delong2015]_ in the NVT ensemble, including both translational and rotational motion.

----

* **type**: ``Brownian``, ``EulerMaruyamaRigidBody``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the implicit solvent :math:`[mass/(distance \cdot time)]` (optional if translational and rotational self-diffusion are defined in particle data)

Example:

.. code-block::

   "eulerMaruyamaRigidBody":{
     "type":["Brownian","EulerMaruyamaRigidBody"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0
     }
   }

.. note::
   If ``viscosity`` is provided, the integrator will compute both translational and rotational self-diffusion coefficients for each particle based on their radii.

.. warning::
   This integrator assumes that particle masses, radii, and orientations (as quaternions) are defined in the particle data.

.. [Delong2015] Delong, Steven, Florencio Balboa Usabiaga, and Aleksandar Donev. "Brownian dynamics of confined rigid bodies." The Journal of chemical physics 143, no. 14 (2015).
