GyrationRadius
--------------

The GyrationRadius step calculates and records the radius of gyration of the system or a specified group of particles.

Output format:

.. code-block::

   # step GyrationRadius
   0 2.345
   100 2.456
   200 2.567

----

* **type**: ``GeometricalMeasure``, ``GyrationRadius``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "gyrationRadius":{
     "type":["GeometricalMeasure","GyrationRadius"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "gyration_radius.dat"
     }
   }

.. note::
   The output file will contain the step number and the radius of gyration value.

.. note::
   A group can be specified by using the ``group`` parameter. The radius of gyration will be calculated for the specified group of particles.

.. tip::
   The radius of gyration is particularly useful for characterizing the size and shape of polymers or other complex molecular structures.
