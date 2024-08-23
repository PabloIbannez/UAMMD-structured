Time
-----

When ``Time`` is selected a set of global variables are defined and used by the other components.
When this component is selected, time is assumed to be the driving variable of the simulation (which is the case for most simulations).


* **type**: ``Fundamental``, ``Time``.
* **parameters**:

  * ``currentStep`` : ``unsigned long long int``, *optional*, default: 0.

  * ``simulationTime`` : ``float``, *optional*, default: 0.0.

  * ``timeStep`` : ``float``, *optional*.

* **data**: ``None``.

The ``currentStep``, ``simulationTime`` and ``timeStep`` are usually not set by the user, but are updated by the simulation engine.

----

Example:

.. code-block:: json

   "entryName": {
     "type": ["Fundamental", "Time"],
     "parameters": {
       "currentStep": 0,
       "timeStep": 0.001,
       "simulationTime": 0.0
     }
   }
