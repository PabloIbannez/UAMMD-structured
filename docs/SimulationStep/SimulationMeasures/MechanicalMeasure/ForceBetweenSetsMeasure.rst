ForceBetweenSetsMeasure
-----------------------

The ForceBetweenSetsMeasure step calculates and records the total force between defined sets of particles.

----

* **type**: ``MechanicalMeasure``, ``ForceBetweenSetsMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``name``: ``string``: Name of the particle set
  * ``id_list``: ``int[]``: List of particle IDs in the set

Example:

.. code-block::

   "forceBetweenSets":{
     "type":["MechanicalMeasure","ForceBetweenSetsMeasure"],
     "parameters":{
       "outputFilePath": "force_between_sets.dat"
     },
     "data":[
       {"name": "set1", "id_list": [1, 2, 3, 4, 5]},
       {"name": "set2", "id_list": [6, 7, 8, 9, 10]}
     ]
   }

.. note::
   The output file will contain the total force between each pair of defined sets.

.. tip::
   This measure is particularly useful for studying interactions between different parts of a system, such as protein domains or nanoparticle clusters.
