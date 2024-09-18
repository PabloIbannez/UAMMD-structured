Quasi2D
-------

The Quasi2D integrator implements Brownian Dynamics with Hydrodynamic Interactions for quasi-two-dimensional systems (periodic in two dimensions, confined in the third). This integrator is a wrapper around the UAMMD BDHI::Quasi2D integrator.

For more details on the underlying method, please refer to the `UAMMD BDHI documentation <https://uammd.readthedocs.io/en/latest/Integrator/BrownianHydrodynamics.html>`_.

----

* **type**: ``BDHIDoublyPeriodic``, ``Quasi2D``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[mass/(distance \cdot time)]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]` (optional, uses particle radius if not specified)
  * ``tolerance``: ``real``: Tolerance for the iterative solver (default: 1e-3)

Example:

.. code-block::

   "quasi2d":{
     "type":["BDHIDoublyPeriodic","Quasi2D"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "hydrodynamicRadius": 0.5,
       "tolerance": 1e-4
     }
   }

.. note::
   This integrator is suitable for systems with periodic boundary conditions in two dimensions and confinement in the third dimension.

.. warning::
   This integrator assumes that particle masses and radii are defined in the particle data.
