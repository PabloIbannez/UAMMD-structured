ForceCouplingMethod
-------------------

The ForceCouplingMethod integrator implements Brownian Dynamics with Hydrodynamic Interactions using the Force Coupling Method for triply periodic systems. This integrator is a wrapper around the UAMMD BDHI::EulerMaruyama<BDHI::FCM> integrator.

For more details on the underlying method, please refer to the `UAMMD BDHI documentation <https://uammd.readthedocs.io/en/latest/Integrator/BrownianHydrodynamics.html>`_.

----

* **type**: ``BDHITriplePeriodic``, ``ForceCouplingMethod``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[mass/(distance \cdot time)]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]` (default: 1.0)
  * ``tolerance``: ``real``: Tolerance for the iterative solver (default: 1e-3)

Example:

.. code-block::

   "fcm":{
     "type":["BDHITriplePeriodic","ForceCouplingMethod"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "hydrodynamicRadius": 0.5,
       "tolerance": 1e-4
     }
   }

.. note::
   This integrator is suitable for systems with periodic boundary conditions in all three dimensions.

.. warning::
   This integrator assumes that particle masses are defined in the particle data.
