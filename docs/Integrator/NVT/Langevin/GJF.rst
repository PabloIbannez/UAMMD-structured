GJF
---

The GJF (Grønbech-Jensen/Farago) integrator implements a Langevin dynamics algorithm for simulating particles in the NVT ensemble. This integrator is a custom implementation in UAMMD-structured.

----

* **type**: ``Langevin``, ``GJF``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``frictionConstant``: ``real``: Friction constant for the Langevin thermostat :math:`[1/time]`
  * ``resetVelocities``: ``bool``: Whether to reset velocities at the start of the simulation (default: true)

Example:

.. code-block::

   "gjf":{
     "type":["Langevin","GJF"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "frictionConstant": 1.0,
       "resetVelocities": true
     }
   }

.. note::
   The GJF integrator is known for its improved accuracy in sampling configurational quantities compared to other Langevin integrators.

.. warning::
   This integrator assumes that the particle mass is defined in the particle data.
