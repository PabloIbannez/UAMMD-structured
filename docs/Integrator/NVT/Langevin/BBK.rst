BBK
---

The BBK (Brünger–Brooks–Karplus) [Brunger1984]_ [Pastor1988]_ [Loncharich1992]_ integrator implements a Langevin dynamics algorithm for simulating particles in the NVT ensemble. This integrator is a custom implementation in UAMMD-structured.

----

* **type**: ``Langevin``, ``BBK``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``frictionConstant``: ``real``: Friction constant for the Langevin thermostat :math:`[1/time]` (optional if defined in particle data)
  * ``stopTransRotSteps``: ``ullint``: Number of steps between removing global translation and rotation (default: 0, disabled)
  * ``resetVelocities``: ``bool``: Whether to reset velocities at the start of the simulation (default: true)

Example:

.. code-block::

   "bbk":{
     "type":["Langevin","BBK"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "frictionConstant": 1.0,
       "stopTransRotSteps": 1000,
       "resetVelocities": true
     }
   }

.. note::
   If ``frictionConstant`` is not provided, the integrator will use the friction constant defined in the particle data.

.. warning::
   This integrator assumes that the particle mass is defined in the particle data.

.. [Brunger1984] Brünger, Axel, Charles L. Brooks III, and Martin Karplus. "Stochastic boundary conditions for molecular dynamics simulations of ST2 water." Chemical physics letters 105, no. 5 (1984): 495-500.
.. [Pastor1988] Pastor, Richard W., Bernard R. Brooks, and Attila Szabo. "An analysis of the accuracy of Langevin and molecular dynamics algorithms." Molecular Physics 65, no. 6 (1988): 1409-1419.
.. [Loncharich1992] Loncharich, Richard J., Bernard R. Brooks, and Richard W. Pastor. "Langevin dynamics of peptides: The frictional dependence of isomerization rates of N‐acetylalanyl‐N′‐methylamide." Biopolymers: Original Research on Biomolecules 32, no. 5 (1992): 523-535.
