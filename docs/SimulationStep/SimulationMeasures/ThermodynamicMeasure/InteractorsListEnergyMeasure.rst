InteractorsListEnergyMeasure
----------------------------

The InteractorsListEnergyMeasure step calculates and records the energy contributions from specified interactors over time.

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
       "outputFilePath": "interactors_energy.dat",
       "interactorsList": ["LennardJones", "Coulomb", "Bonds"]
     }
   }

.. note::
   The output file will contain columns for the step number and the energy contribution from each specified interactor.

.. tip::
   This measure is useful for analyzing the relative contributions of different interaction types to the total system energy.
