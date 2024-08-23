AnglesMeasure
-------------

The AnglesMeasure step calculates and records the angles formed by specified triplets of particles over time.

----

* **type**: ``ParticlesListMeasure``, ``AnglesMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``id_i``: ``int``: ID of the first particle in the angle
  * ``id_j``: ``int``: ID of the central particle in the angle
  * ``id_k``: ``int``: ID of the third particle in the angle

Example:

.. code-block::

   "anglesMeasure":{
     "type":["ParticlesListMeasure","AnglesMeasure"],
     "parameters":{
       "outputFilePath": "angles.dat"
     },
     "labels":["id_i", "id_j", "id_k"],
     "data":[
       [1, 2, 3],
       [2, 3, 4],
       [3, 4, 5]
     ]
   }

.. note::
   The output file will contain the step number and the angle for each specified triplet of particles.

.. tip::
   This measure is useful for analyzing bond angles, molecular conformations, or the relative orientation of specific groups of particles in the system.
