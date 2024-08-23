Minimization
============

SteepestDescent
---------------

The SteepestDescent integrator implements a simple energy minimization algorithm using the method of steepest descent.

----

* **type**: ``Minimization``, ``SteepestDescent``
* **parameters**:

  * ``h``: ``real``: Step size for the minimization
  * ``maxObjectiveForce``: ``real``: Maximum allowed force to consider the system minimized
  * ``nStepsPrintProgress``: ``ullint``: Number of steps between progress reports (default: 0, disabled)

Example:

.. code-block::

   "steepestDescent":{
     "type":["Minimization","SteepestDescent"],
     "parameters":{
       "h": 0.01,
       "maxObjectiveForce": 1e-4,
       "nStepsPrintProgress": 1000
     }
   }

.. note::
   The minimization stops when the maximum force in the system becomes smaller than ``maxObjectiveForce``.

.. warning::
   This integrator modifies particle positions but does not update velocities.
