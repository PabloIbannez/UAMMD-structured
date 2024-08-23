LambdaActivation
----------------

The LambdaActivation step controls the lambda parameter for alchemical transformations or free energy calculations. It allows for gradual changes in lambda according to a predefined schedule.

----

* **type**: ``FlowControl``, ``LambdaActivation``
* **parameters**:

  * ``lambdaValueStep``: ``int``: Number of steps between lambda value changes
  * ``lambdaValues``: ``real[]``: List of lambda values to cycle through

Example:

.. code-block::

   "lambdaActivation":{
     "type":["FlowControl","LambdaActivation"],
     "parameters":{
       "lambdaValueStep": 1000,
       "lambdaValues": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
     }
   }

.. note::
   This step is useful for free energy calculations and alchemical transformations.

.. warning::
   Ensure that your forcefield and integrator support lambda-dependent interactions when using this step.
