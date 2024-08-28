ContactsMeasure
---------------

The ContactsMeasure step calculates and records the distances between specified pairs of particles over time.

Output format:

.. code-block::

   # Step 1
   id1 id2 distance
   id3 id4 distance
   ...
   # Step 2
   id1 id2 distance
   id3 id4 distance
   ...

----

* **type**: ``ParticlesListMeasure``, ``ContactsMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``id_i``: ``int``: ID of the first particle in the pair
  * ``id_j``: ``int``: ID of the second particle in the pair

Example:

.. code-block::

   "contactsMeasure":{
     "type":["ParticlesListMeasure","ContactsMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "contacts.dat"
     },
     "labels":["id_i", "id_j"],
     "data":[
       [1, 5],
       [2, 6],
       [3, 7]
     ]
   }

.. note::
   The output file will contain the step number and the distance for each specified pair of particles.

.. tip::
   This measure is useful for tracking specific interactions or distances in the system, such as monitoring protein contacts or ligand-receptor distances.
