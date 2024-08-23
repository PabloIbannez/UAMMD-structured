HessianMeasure
--------------

The HessianMeasure step calculates and records the Hessian matrix of the system, providing information about the curvature of the potential energy surface.

----

* **type**: ``MechanicalMeasure``, ``HessianMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``mode``: ``string``: Calculation mode, either "Analytical" or "Numerical" (default: "Analytical")
  * ``outputPrecision``: ``int``: Number of decimal places in the output (default: 10)

Example:

.. code-block::

   "hessianMeasure":{
     "type":["MechanicalMeasure","HessianMeasure"],
     "parameters":{
       "outputFilePath": "hessian_measure.dat",
       "mode": "Analytical",
       "outputPrecision": 12
     }
   }

.. note::
   The output file will contain the full Hessian matrix elements for each pair of particles.

.. warning::
   The Hessian calculation can be computationally expensive, especially for large systems.

.. tip::
   The Hessian matrix is useful for normal mode analysis and studying system stability.
