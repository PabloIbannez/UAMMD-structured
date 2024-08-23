Verlet
------

The Verlet integrator for SPH (Smoothed Particle Hydrodynamics) implements a velocity Verlet algorithm combined with SPH interactions. This integrator is a wrapper around the UAMMD VerletNVE integrator combined with an SPH potential.

For more details on the underlying method, please refer to the `UAMMD VerletNVE documentation <https://uammd.readthedocs.io/en/latest/Integrators.html#verletnve>`_ and `UAMMD SPH documentation <https://uammd.readthedocs.io/en/latest/Interactors.html#sph>`_.

----

* **type**: ``SPH``, ``Verlet``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Initial temperature for velocity initialization :math:`[energy]`
  * ``support``: ``real``: Support radius for SPH kernel :math:`[distance]`
  * ``viscosity``: ``real``: Kinematic viscosity of the fluid :math:`[distance^2/time]`
  * ``gasStiffness``: ``real``: Gas stiffness constant :math:`[energy]`
  * ``restDensity``: ``real``: Rest density of the fluid :math:`[mass/distance^3]`

Example:

.. code-block::

   "sphVerlet":{
     "type":["SPH","Verlet"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "support": 1.0,
       "viscosity": 0.1,
       "gasStiffness": 1.0,
       "restDensity": 1.0
     }
   }

.. note::
   This integrator automatically initializes particle velocities according to the specified temperature if they are not already set.

.. warning::
   This integrator assumes that particle masses are defined in the particle data. The SPH formulation also requires particle densities to be computed and stored.
