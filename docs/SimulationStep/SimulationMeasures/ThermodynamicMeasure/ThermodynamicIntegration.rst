ThermodynamicIntegration
------------------------

The ThermodynamicIntegration step performs thermodynamic integration calculations, useful for free energy calculations and alchemical transformations.

The output file contains the value of:

.. math::

   \frac{\partial U}{\partial \lambda}

where U is the potential energy of the system and λ is the lambda value, the format is:

.. code-block::

   # 0.0
   -500.5
   -501.2
   -499.8
   ...
   # 0.1
   -450.3
   -451.1
   -449.7
   ...

where the first line is the lambda value and the following lines are the value of the derivative of the potential energy with respect to lambda.

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
       "intervalStep": 10000,
       "outputFilePath": "thermo_integration.dat",
       "intervalStep": 10000,
       "stepLambda": 10000,
       "lambdaValues": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
     }
   }

.. note::
   The output file will contain the lambda value and the corresponding dU/dλ for each step.

.. warning::
   Ensure that your forcefield and integrator support lambda-dependent interactions when using this measure.

.. tip::
   Thermodynamic integration is a powerful method for calculating free energy differences between two states of a system.
