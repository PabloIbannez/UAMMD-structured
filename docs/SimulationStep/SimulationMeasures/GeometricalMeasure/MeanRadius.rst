MeanRadius
----------

The MeanRadius step calculates and records the mean radius of the system, defined as the average distance of particles from the center of mass.

----

* **type**: ``GeometricalMeasure``, ``MeanRadius``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "meanRadius":{
     "type":["GeometricalMeasure","MeanRadius"],
     "parameters":{
       "outputFilePath": "mean_radius.dat"
     }
   }

.. note::
   The output file will contain the step number and the mean radius value.

.. tip::
   This measure can be useful for monitoring the overall size or compactness of the system over time.
