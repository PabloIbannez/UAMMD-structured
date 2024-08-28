DistancesMeasure
----------------

The DistancesMeasure step calculates and records the distances between specified pairs of particles over time.

Output format:

.. code-block::

   # Step 1
   distance1 distance2 distance3 ...
   # Step 2
   distance1 distance2 distance3 ...
   # Step 3
   distance1 distance2 distance3 ...
   ...

----

* **type**: ``ParticlesListMeasure``, ``DistancesMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``id_i``: ``int``: ID of the first particle in the pair
  * ``id_j``: ``int``: ID of the second particle in the pair

Example:

.. code-block::

   "distancesMeasure":{
     "type":["ParticlesListMeasure","DistancesMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "distances.dat"
     },
     "labels":["id_i", "id_j"],
     "data":[
       [1, 5],
       [2, 6],
       [3, 7]
     ]
   }

.. note::
   The output file will contain the distances for all specified pairs of particles at each step.

.. tip::
   This measure is similar to ContactsMeasure but focuses specifically on distances. It's useful for tracking bond lengths, intermolecular distances, or monitoring separation between specific particles.
