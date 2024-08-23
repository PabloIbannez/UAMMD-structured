StressMeasure
-------------

The StressMeasure step calculates and records the stress tensor for each particle in the system, providing detailed information about local stresses.

----

* **type**: ``MechanicalMeasure``, ``StressMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``radiusCutOff``: ``real``: Cutoff radius for stress calculations :math:`[distance]`

Example:

.. code-block::

   "stressMeasure":{
     "type":["MechanicalMeasure","StressMeasure"],
     "parameters":{
       "outputFilePath": "stress_measure.dat",
       "radiusCutOff": 2.5
     }
   }

.. note::
   The output file will contain particle positions, radii, types, volumes, and the full stress tensor for each particle.

.. tip::
   This measure is particularly useful for studying local stress distributions and identifying stress concentrations in the system.
