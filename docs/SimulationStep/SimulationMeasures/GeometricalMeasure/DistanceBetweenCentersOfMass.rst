DistanceBetweenCentersOfMass
----------------------------

The DistanceBetweenCentersOfMass step calculates and records the distances between the centers of mass of specified groups of particles.

Output format:

.. code-block::

   # step distance_group_pair_1 distance_group_pair_2 distance_group_pair_3 ...
   0 5.678 10.123 15.456
   100 5.789 10.234 15.567
   200 5.890 10.345 15.678

----

* **type**: ``GeometricalMeasure``, ``DistanceBetweenCentersOfMass``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``name``: ``string``: Name of the group
  * ``id_list``: ``int[]``: List of particle IDs in the group

Example:

.. code-block::

   "comDistances":{
     "type":["GeometricalMeasure","DistanceBetweenCentersOfMass"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "com_distances.dat"
     },
     "labels":["idSet_i","idSet_j"],
     "data":[
       [[1, 2, 3, 4, 5],[6, 7, 8, 9, 10]],
       [[11, 12, 13, 14, 15],[16, 17, 18, 19, 20]],
       [[21, 22, 23, 24, 25],[26, 27, 28, 29, 30]]
     ]
   }

.. note::
   The output file will contain the step number and the distances between all pairs of specified groups.

.. tip::
   This measure is useful for tracking the relative motion of different parts of the system, such as in protein-ligand interactions or polymer aggregation studies.
