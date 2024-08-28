PairwiseForceMeasure
--------------------

The PairwiseForceMeasure step calculates and records the pairwise forces between particles or the total force on each particle.

Output format:

.. code-block::

   # For mode "Pairwise_force":
   #From     To       Fx                  Fy                  Fz
   1         2        1.2345678901        -0.9876543210       0.5432109876
   1         3        -0.1111111111       0.2222222222        -0.3333333333
   ...

   # For mode "Total_force":
   #Id       Fx                  Fy                  Fz
   1         1.1234567890        -0.7654321098       0.2098765432
   2         -0.9876543210       0.8765432109        -0.7654321098
   ...

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
       "intervalStep": 10000,
       "outputFilePath": "pairwise_force.dat",
       "mode": "Pairwise_force"
     }
   }

.. note::
   In "Pairwise_force" mode, the output file will contain forces between each pair of interacting particles. In "Total_force" mode, it will contain the total force on each particle.

.. tip::
   This measure is useful for analyzing force distributions and identifying strong interactions in the system.
