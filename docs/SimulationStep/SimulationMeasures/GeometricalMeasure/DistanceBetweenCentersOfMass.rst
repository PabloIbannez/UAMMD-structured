DistanceBetweenCentersOfMass
----------------------------

The DistanceBetweenCentersOfMass step calculates and records the distances between the centers of mass of specified groups of particles.

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
       "outputFilePath": "com_distances.dat"
     },
     "data":[
       {"name": "group1", "id_list": [1, 2, 3, 4, 5]},
       {"name": "group2", "id_list": [6, 7, 8, 9, 10]}
     ]
   }

.. note::
   The output file will contain the step number and the distances between all pairs of specified groups.

.. tip::
   This measure is useful for tracking the relative motion of different parts of the system, such as in protein-ligand interactions or polymer aggregation studies.
