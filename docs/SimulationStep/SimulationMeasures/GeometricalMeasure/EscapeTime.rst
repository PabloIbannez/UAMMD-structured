EscapeTime
----------

The EscapeTime step records the time at which particles escape from a defined region in the simulation box.

Output format:

.. code-block::

   0 3
   0 7
   100 12
   200 18
   300 25

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
       "intervalStep": 10000,
       "outputFilePath": "escape_times.dat"
     },
     "labels":["normalVector","independentVector"],
     "data":[
        [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],[[5, 0, 0], [0, 5, 0], [0, 0, 5]]],
        [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],[[2,1,0], [0, 2, 0], [0, 0, 2]]]
     ]
   }

.. note::
   The output file will contain the escape time and particle ID for each particle that leaves the defined region.

.. tip::
   This measure is useful for studying confinement effects, particle trapping, or escape rates in various systems.
