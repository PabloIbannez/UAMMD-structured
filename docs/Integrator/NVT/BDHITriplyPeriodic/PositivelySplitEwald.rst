PositivelySplitEwald
--------------------

The PositivelySplitEwald integrator implements Brownian Dynamics with Hydrodynamic Interactions using the Positively Split Ewald method for triply periodic systems. This integrator is a wrapper around the UAMMD BDHI::EulerMaruyama<BDHI::PSE> integrator.

For more details on the underlying method, please refer to the `UAMMD BDHI documentation <https://uammd.readthedocs.io/en/latest/Integrator/BrownianHydrodynamics.html>`_.

----

* **type**: ``BDHITriplePeriodic``, ``PositivelySplitEwald``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[mass/(distance \cdot time)]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]` (default: 1.0)
  * ``tolerance``: ``real``: Tolerance for the Ewald sum (default: 1e-3)
  * ``psi``: ``real``: Ewald splitting parameter (default: 0.5)

Example:

.. code-block::

   "pse":{
     "type":["BDHITriplePeriodic","PositivelySplitEwald"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "hydrodynamicRadius": 0.5,
       "tolerance": 1e-4,
       "psi": 0.6
     }
   }

.. note::
   This integrator is suitable for systems with periodic boundary conditions in all three dimensions.

.. warning::
   This integrator assumes that particle masses are defined in the particle data.
