LambdaCycle
-----------

The LambdaCycle step implements a cyclic lambda schedule for alchemical transformations or free energy calculations. It allows for repeated cycles of lambda changes, useful for measuring free energy of particles immersed in a solvent. The cycle consists of a series of steps for each lambda value, followed by a pause step before the next cycle begins. The cycle is repeated until the simulation is complete.

.. code-block::

    lambda=1                        ...................                               ...................
    lambdaValues[n]              ...                                               ...
        .......               ...                                               ...
        .......            ...                                               ...
    lambdaValues[0]     ...                                               ...
    lambda=0                                           ...................                               ...................
    startStep ---------><-><-><-><-><---measureStep---><----pauseStep----><-><-><-><-><---measureStep---><----pauseStep---->
                         |
                         |
                         L-> activationStep

----

* **type**: ``FlowControl``, ``LambdaCycle``
* **parameters**:

  * ``activationStep``: ``int``: Number of steps for each lambda activation phase
  * ``measureStep``: ``int``: Number of steps for measurement at each lambda value
  * ``pauseStep``: ``int``: Number of steps for pause between cycles
  * ``lambdaValues``: ``real[]``: List of lambda values to cycle through

Example:

.. code-block::

   "lambdaCycle":{
     "type":["FlowControl","LambdaCycle"],
     "parameters":{
       "activationStep": 1000,
       "measureStep": 5000,
       "pauseStep": 2000,
       "lambdaValues": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
     }
   }

.. note::
   This step is particularly useful for enhanced sampling techniques like lambda dynamics.

.. warning::
   Careful consideration of the cycle parameters is necessary to ensure proper equilibration and sampling at each lambda value.
