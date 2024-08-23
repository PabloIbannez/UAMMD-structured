MeasureMeanMagnetization
------------------------

The MeasureMeanMagnetization step calculates and records the mean magnetization of the system over time, normalized by the maximum possible magnetization.

----

* **type**: ``MagneticMeasure``, ``MeasureMeanMagnetization``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``startStep``: ``int``: Step at which to start the measurement (default: 0)

Example:

.. code-block::

   "meanMagnetization":{
     "type":["MagneticMeasure","MeasureMeanMagnetization"],
     "parameters":{
       "outputFilePath": "mean_magnetization.dat",
       "startStep": 1000
     }
   }

.. note::
   The output file will contain the time and the x, y, z components of the normalized mean magnetization vector.

.. tip::
   This measure provides insight into the degree of magnetic alignment in the system, with values ranging from 0 (random orientation) to 1 (perfect alignment).
