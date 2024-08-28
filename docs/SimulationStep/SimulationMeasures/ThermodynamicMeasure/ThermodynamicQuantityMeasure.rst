ThermodynamicQuantityMeasure
----------------------------

The ThermodynamicQuantityMeasure step calculates and records various thermodynamic quantities of the system over time, including energy components, temperature, and virial.

Output format:

.. code-block::

   #Step ParticleNumber Volume Energy(Interactor1) Energy(Interactor2) ... KineticEnergy TotalPotentialEnergy TotalEnergy Temperature Virial
   0     1000           1000.0 -500.5              -200.3              ... 300.2         -700.8               -400.6      298.15      -150.3
   100   1000           1000.0 -502.1              -201.5              ... 301.8         -703.6               -401.8      298.20      -151.2
   ...

----

* **type**: ``ThermodynamicMeasure``, ``ThermodynamicQuantityMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "thermodynamicQuantities":{
     "type":["ThermodynamicMeasure","ThermodynamicQuantityMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "thermo_quantities.dat"
     }
   }

.. note::
   The output file will contain columns for step number, particle number, volume, energy components (for each interactor), kinetic energy, total potential energy, total energy, temperature, and virial.

.. tip::
   This measure provides a comprehensive overview of the system's thermodynamic state and is useful for monitoring equilibration and system stability.
