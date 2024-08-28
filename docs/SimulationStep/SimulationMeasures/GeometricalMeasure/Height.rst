Height
------

The Height step calculates and records the average height of the highest particles in the system.

----

* **type**: ``GeometricalMeasure``, ``Height``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``particleNumberAverage``: ``int``: Number of highest particles to average (default: 1)

Example:

.. code-block::

   "systemHeight":{
     "type":["GeometricalMeasure","Height"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "system_height.dat",
       "particleNumberAverage": 10
     }
   }

.. note::
   The output file will contain the step number and the average height of the specified number of highest particles.

.. tip::
   This measure can be useful for studying surface phenomena or monitoring the expansion/compression of the system in the vertical direction.
