StressMeasure
-------------

The StressMeasure step calculates and records the stress tensor for each particle in the system, providing detailed information about local stresses. The employed procedure is described in detail in [kononova2018]_.

Output format:

.. code-block::

   #Lx=10.0;Ly=10.0;Lz=10.0;
   x y z          radius type volume σxx   σxy   σxz   σyx   σyy   σyz   σzx   σzy   σzz
   1.23 4.56 7.89 0.5    1    0.5236 10.0  0.1   0.2   0.1   20.0  0.3   0.2   0.3   30.0
   2.34 5.67 8.90 0.6    2    0.9048 15.0  0.2   0.3   0.2   25.0  0.4   0.3   0.4   35.0
   ...

----

* **type**: ``MechanicalMeasure``, ``StressMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``radiusCutOff``: ``real``: Cutoff radius for stress calculations :math:`[distance]`

Example:

.. code-block::

   "stressMeasure":{
     "type":["MechanicalMeasure","StressMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "stress_measure.dat",
       "radiusCutOff": 2.5
     }
   }

.. note::
   The output file will contain particle positions, radii, types, volumes, and the full stress tensor for each particle.

.. tip::
   This measure is particularly useful for studying local stress distributions and identifying stress concentrations in the system.

.. [kononova2018] Kononova, Olga, Farkhad Maksudov, Kenneth A. Marx, and Valeri Barsegov. "TensorCalculator: exploring the evolution of mechanical stress in the CCMV capsid." Journal of Physics: Condensed Matter 30, no. 4 (2018): 044006.
