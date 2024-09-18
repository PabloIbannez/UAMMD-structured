Cholesky
--------

The Cholesky integrator implements Brownian Dynamics with Hydrodynamic Interactions for open boundary systems using the Cholesky decomposition method. This integrator is a wrapper around the UAMMD BDHI::EulerMaruyama<BDHI::Cholesky> integrator.

For more details on the underlying method, please refer to the `UAMMD BDHI documentation <https://uammd.readthedocs.io/en/latest/Integrator/BrownianHydrodynamics.html>`_.

----

* **type**: ``BDHIOpenBoundary``, ``Cholesky``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[mass/(distance \cdot time)]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]` (optional, uses particle radius if not specified)

Example:

.. code-block::

   "cholesky":{
     "type":["BDHIOpenBoundary","Cholesky"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "hydrodynamicRadius": 0.5
     }
   }

.. note::
   This integrator is suitable for systems with open boundary conditions in all directions. It provides exact hydrodynamic interactions but scales poorly with system size.

.. warning::
   This integrator assumes that particle masses and radii are defined in the particle data. Due to its computational complexity, it is only suitable for small systems (typically less than 1000 particles).
