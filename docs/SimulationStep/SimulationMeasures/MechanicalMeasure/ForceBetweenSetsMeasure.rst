ForceBetweenSetsMeasure
-----------------------

The ForceBetweenSetsMeasure step calculates and records the total force between defined sets of particles.

Output format:

.. code-block::

   #
   From           To             Fx                  Fy                  Fz
   set1           set2           -1.2345678901       0.9876543210        0.5432109876
   set2           set1           1.2345678901        -0.9876543210       -0.5432109876
   set1           set3           0.1111111111        -0.2222222222       0.3333333333
   set3           set1           -0.1111111111       0.2222222222        -0.3333333333
   ...

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
       "intervalStep": 10000,
       "outputFilePath": "force_between_sets.dat"
     },
     "labels":["name","id_list"],
     "data":[
       ["set1", [1, 2, 3, 4, 5]],
       ["set2", [6, 7, 8, 9, 10]]
     ]
   }

.. note::
   The output file will contain the total force between each pair of defined sets.

.. tip::
   This measure is particularly useful for studying interactions between different parts of a system, such as protein domains or nanoparticle clusters.
