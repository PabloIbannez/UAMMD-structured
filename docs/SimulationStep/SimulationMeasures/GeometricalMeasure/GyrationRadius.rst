GyrationRadius
--------------

The GyrationRadius step calculates and records the radius of gyration of the system or a specified group of particles.

----

* **type**: ``GeometricalMeasure``, ``GyrationRadius``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "gyrationRadius":{
     "type":["GeometricalMeasure","GyrationRadius"],
     "parameters":{
       "outputFilePath": "gyration_radius.dat"
     }
   }

.. note::
   The output file will contain the step number and the radius of gyration value.

.. tip::
   The radius of gyration is particularly useful for characterizing the size and shape of polymers or other complex molecular structures.
