CenterOfMassPosition
--------------------

The CenterOfMassPosition step calculates and records the center of mass position of the system or a specified group of particles over time.

----

* **type**: ``GeometricalMeasure``, ``CenterOfMassPosition``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "comPosition":{
     "type":["GeometricalMeasure","CenterOfMassPosition"],
     "parameters":{
       "outputFilePath": "com_position.dat"
     }
   }

.. note::
   The output file will contain the step number and the x, y, z coordinates of the center of mass.

.. tip::
   This measure is useful for tracking the overall motion of the system or specific groups of particles.
