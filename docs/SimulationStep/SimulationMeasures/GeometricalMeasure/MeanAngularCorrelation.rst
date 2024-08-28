MeanAngularCorrelation
----------------------

The MeanAngularCorrelation step calculates and records the mean angular correlation of particles' orientations relative to their initial orientations.

Output format:

.. code-block::

   # step MAC
   0 1.000
   100 0.987
   200 0.954

----

* **type**: ``GeometricalMeasure``, ``MeanAngularCorrelation``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "angularCorrelation":{
     "type":["GeometricalMeasure","MeanAngularCorrelation"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "angular_correlation.dat"
     }
   }

.. note::
   The output file will contain the step number and the mean angular correlation value.

.. tip::
   This measure is particularly useful for studying the rotational dynamics of particles in the system.
