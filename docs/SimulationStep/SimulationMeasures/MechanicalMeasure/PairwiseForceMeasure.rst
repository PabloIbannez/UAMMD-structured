PairwiseForceMeasure
--------------------

The PairwiseForceMeasure step calculates and records the pairwise forces between particles or the total force on each particle.

----

* **type**: ``MechanicalMeasure``, ``PairwiseForceMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``mode``: ``string``: Calculation mode, either "Pairwise_force" or "Total_force" (default: "Pairwise_force")

Example:

.. code-block::

   "pairwiseForceMeasure":{
     "type":["MechanicalMeasure","PairwiseForceMeasure"],
     "parameters":{
       "outputFilePath": "pairwise_force.dat",
       "mode": "Pairwise_force"
     }
   }

.. note::
   In "Pairwise_force" mode, the output file will contain forces between each pair of interacting particles. In "Total_force" mode, it will contain the total force on each particle.

.. tip::
   This measure is useful for analyzing force distributions and identifying strong interactions in the system.
