EscapeTime
----------

The EscapeTime step records the time at which particles escape from a defined region in the simulation box.

----

* **type**: ``GeometricalMeasure``, ``EscapeTime``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``normalVector``: ``real3[]``: List of normal vectors defining the escape planes
  * ``independentVector``: ``real3[]``: List of points on the escape planes

Example:

.. code-block::

   "escapeTime":{
     "type":["GeometricalMeasure","EscapeTime"],
     "parameters":{
       "outputFilePath": "escape_times.dat"
     },
     "data":{
       "normalVector": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
       "independentVector": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
     }
   }

.. note::
   The output file will contain the escape time and particle ID for each particle that leaves the defined region.

.. tip::
   This measure is useful for studying confinement effects, particle trapping, or escape rates in various systems.
