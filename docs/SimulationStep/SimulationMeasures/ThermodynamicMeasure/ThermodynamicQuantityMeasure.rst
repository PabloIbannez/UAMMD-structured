ThermodynamicQuantityMeasure
----------------------------

The ThermodynamicQuantityMeasure step calculates and records various thermodynamic quantities of the system over time, including energy components, temperature, and virial.

----

* **type**: ``ThermodynamicMeasure``, ``ThermodynamicQuantityMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "thermodynamicQuantities":{
     "type":["ThermodynamicMeasure","ThermodynamicQuantityMeasure"],
     "parameters":{
       "outputFilePath": "thermo_quantities.dat"
     }
   }

.. note::
   The output file will contain columns for step number, particle number, volume, energy components (for each interactor), kinetic energy, total potential energy, total energy, temperature, and virial.

.. tip::
   This measure provides a comprehensive overview of the system's thermodynamic state and is useful for monitoring equilibration and system stability.
