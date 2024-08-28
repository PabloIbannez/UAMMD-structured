InteractorsListEnergyMeasure
----------------------------

The InteractorsListEnergyMeasure step calculates and records the energy contributions from specified interactors over time.

Output format:

.. code-block::

   # Step     Interactor1 Interactor2 Interactor3 ...
   0          -500.5      -200.3      -100.2      ...
   100        -502.1      -201.5      -101.8      ...
   200        -501.8      -200.9      -100.5      ...
   ...

----

* **type**: ``ThermodynamicMeasure``, ``InteractorsListEnergyMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``interactorsList``: ``string[]``: List of interactor names to measure (optional, if not provided, all interactors will be measured)

Example:

.. code-block::

   "interactorsEnergy":{
     "type":["ThermodynamicMeasure","InteractorsListEnergyMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "interactors_energy.dat",
       "interactorsList": ["LennardJones", "Coulomb", "Bonds"]
     }
   }

.. note::
   The output file will contain columns for the step number and the energy contribution from each specified interactor.

.. tip::
   This measure is useful for analyzing the relative contributions of different interaction types to the total system energy.
