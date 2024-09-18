Verlet
------

The Verlet integrator for DPD (Dissipative Particle Dynamics) implements a velocity Verlet algorithm with DPD thermostat. This integrator is a wrapper around the UAMMD VerletNVE integrator combined with a DPD potential.

For more details on the underlying method, please refer to the `UAMMD VerletNVE documentation <https://uammd.readthedocs.io/en/latest/Integrator/MolecularDynamics.html#verletnve>`_ and `UAMMD DPD Potential documentation <https://uammd.readthedocs.io/en/latest/Interactor/PairForces.html#dpd>`_.

----

* **type**: ``DPD``, ``Verlet``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``cutOff``: ``real``: Cutoff distance for DPD interactions :math:`[distance]`
  * ``gamma``: ``real``: Friction coefficient for DPD interactions :math:`[mass/time]`
  * ``alpha``: ``real``: Strength of the conservative force in DPD :math:`[energy/distance]`

Example:

.. code-block::

   "dpdVerlet":{
     "type":["DPD","Verlet"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "cutOff": 1.0,
       "gamma": 4.5,
       "alpha": 25.0
     }
   }

.. note::
   This integrator automatically initializes particle velocities according to the specified temperature if they are not already set.

.. warning::
   This integrator assumes that particle masses are defined in the particle data.
