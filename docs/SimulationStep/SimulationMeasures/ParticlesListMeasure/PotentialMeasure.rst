PotentialMeasure
----------------

The PotentialMeasure step calculates and records the potential energy, forces, and torques for specified particles, broken down by interaction type.

Output format:

.. code-block::

   # id            InteractorAEnergy InteractorAForceX InteractorAForceY InteractorAForceZ InteractorAForceMod InteractorATorqueX InteractorATorqueY InteractorATorqueZ InteractorATorqueMod InteractorBEnergy ...
   # Step 1
   1               energy_A_1        forceX_A_1        forceY_A_1        forceZ_A_1        forceMod_A_1        torqueX_A_1        torqueY_A_1        torqueZ_A_1        torqueMod_A_1        energy_B_1        ...
   2               energy_A_2        forceX_A_2        forceY_A_2        forceZ_A_2        forceMod_A_2        torqueX_A_2        torqueY_A_2        torqueZ_A_2        torqueMod_A_2        energy_B_2        ...
   ...
   # Step 2
   1               energy_A_1        forceX_A_1        forceY_A_1        forceZ_A_1        forceMod_A_1        torqueX_A_1        torqueY_A_1        torqueZ_A_1        torqueMod_A_1        energy_B_1        ...
   2               energy_A_2        forceX_A_2        forceY_A_2        forceZ_A_2        forceMod_A_2        torqueX_A_2        torqueY_A_2        torqueZ_A_2        torqueMod_A_2        energy_B_2        ...
   ...

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
       "intervalStep": 10000,
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
