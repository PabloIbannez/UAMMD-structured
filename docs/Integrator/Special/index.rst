Special
=======

None
----

The None integrator is a special integrator that does not perform any integration steps. It can be used for testing or in situations where you want to manually control the simulation progression.

----

* **type**: ``Special``, ``None``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]` (used only for updating the simulation time)

Example:

.. code-block::

   "noneIntegrator":{
     "type":["Special","None"],
     "parameters":{
       "timeStep": 0.01
     }
   }

.. note::
   This integrator only updates the current step and simulation time without modifying particle positions or velocities.
