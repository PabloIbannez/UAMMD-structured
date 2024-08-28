CenterOfMassPosition
--------------------

The CenterOfMassPosition step calculates and records the center of mass position of the system or a specified group of particles over time.

Output format:

.. code-block::

   # step CenterOfMass_X CenterOffMass_Y CenterOffMass_Z
   0 1.234 5.678 9.012
   100 1.345 5.789 9.123
   200 1.456 5.890 9.234

----

* **type**: ``GeometricalMeasure``, ``CenterOfMassPosition``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "comPosition":{
     "type":["GeometricalMeasure","CenterOfMassPosition"],
     "parameters":{
       "intervalStep": 10000,
       "group": "selected_group",
       "outputFilePath": "com_position.dat"
     }
   }

.. note::
   The output file will contain the step number and the x, y, z coordinates of the center of mass.

.. tip::
   This measure is useful for tracking the overall motion of the system or specific groups of particles.
