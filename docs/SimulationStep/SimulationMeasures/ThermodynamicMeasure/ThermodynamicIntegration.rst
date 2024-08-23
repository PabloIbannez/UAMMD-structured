ThermodynamicIntegration
------------------------

The ThermodynamicIntegration step performs thermodynamic integration calculations, useful for free energy calculations and alchemical transformations.

----

* **type**: ``ThermodynamicMeasure``, ``ThermodynamicIntegration``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``stepLambda``: ``int``: Number of steps between lambda value changes
  * ``lambdaValues``: ``real[]``: List of lambda values to use

Example:

.. code-block::

   "thermodynamicIntegration":{
     "type":["ThermodynamicMeasure","ThermodynamicIntegration"],
     "parameters":{
       "outputFilePath": "thermo_integration.dat",
       "stepLambda": 1000,
       "lambdaValues": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
     }
   }

.. note::
   The output file will contain the lambda value and the corresponding dU/dλ for each step.

.. warning::
   Ensure that your forcefield and integrator support lambda-dependent interactions when using this measure.

.. tip::
   Thermodynamic integration is a powerful method for calculating free energy differences between two states of a system.
