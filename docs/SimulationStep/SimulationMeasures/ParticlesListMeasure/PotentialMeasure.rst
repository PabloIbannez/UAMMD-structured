PotentialMeasure
----------------

The PotentialMeasure step calculates and records the potential energy, forces, and torques for specified particles, broken down by interaction type.

----

* **type**: ``ParticlesListMeasure``, ``PotentialMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``id``: ``int``: ID of the particle to measure

Example:

.. code-block::

   "potentialMeasure":{
     "type":["ParticlesListMeasure","PotentialMeasure"],
     "parameters":{
       "outputFilePath": "potential_measure.dat"
     },
     "labels":["id"],
     "data":[
       [1],
       [2],
       [3],
       [4],
       [5]
     ]
   }

.. note::
   The output file will contain detailed energy, force, and torque information for each specified particle, broken down by interaction type.

.. tip::
   This measure provides in-depth information about the energetics and forces acting on specific particles of interest in the system.
