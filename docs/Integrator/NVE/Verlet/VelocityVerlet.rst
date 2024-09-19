VelocityVerlet
--------------

The VelocityVerlet integrator implements the velocity Verlet algorithm for simulating particles in the NVE ensemble. This integrator is a wrapper around the UAMMD VerletNVE integrator.

For more details on the underlying method, please refer to the `UAMMD VerletNVE documentation <https://uammd.readthedocs.io/en/latest/Integrator/MolecularDynamics.html#molecular-dynamics>`_.

----

* **type**: ``Verlet``, ``VelocityVerlet``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Initial temperature for velocity initialization :math:`[temperature]`

Example:

.. code-block::

   "velocityVerlet":{
     "type":["Verlet","VelocityVerlet"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0
     }
   }

.. note::
   The temperature parameter is only used for initializing velocities at the start of the simulation if they are not already set.
