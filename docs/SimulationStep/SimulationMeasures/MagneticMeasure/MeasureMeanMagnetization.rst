MeasureTotalMagnetization
-------------------------

The MeasureTotalMagnetization step calculates and records the total magnetization of the system over time.

----

* **type**: ``MagneticMeasure``, ``MeasureTotalMagnetization``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``startStep``: ``int``: Step at which to start the measurement (default: 0)

Example:

.. code-block::

   "totalMagnetization":{
     "type":["MagneticMeasure","MeasureTotalMagnetization"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "total_magnetization.dat",
       "startStep": 1000
     }
   }

.. note::
   The output file will contain the time and the x, y, z components of the total magnetization vector.

.. tip::
   This measure is useful for studying the overall magnetic behavior of the system, especially in response to external fields or temperature changes.
