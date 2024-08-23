EulerMaruyama
-------------

The EulerMaruyama integrator implements the Euler-Maruyama [Desmond2001]_ method for simulating Brownian dynamics in the NVT ensemble.

----

* **type**: ``Brownian``, ``EulerMaruyama``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``viscosity``: ``real``: Viscosity of the implicit solvent :math:`[mass/(distance \cdot time)]` (optional if translational self-diffusion is defined in particle data)

Example:

.. code-block::

   "eulerMaruyama":{
     "type":["Brownian","EulerMaruyama"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0
     }
   }

.. note::
   If ``viscosity`` is provided, the integrator will compute the translational self-diffusion coefficients for each particle based on their radii.

.. warning::
   This integrator assumes that particle masses and radii are defined in the particle data.

.. [Desmond2001] Desmond J. Higham. An algorithmic introduction to numerical simulation of stochastic differential equations. SIAM review, 43(3):525–546, 2001.
