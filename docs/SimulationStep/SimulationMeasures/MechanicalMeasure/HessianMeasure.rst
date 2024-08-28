HessianMeasure
--------------

The HessianMeasure step calculates and records the Hessian matrix of the system, providing information about the curvature of the potential energy surface.

Output format:

.. code-block::

   #Id_i     Id_j     Hxx                 Hxy                 Hxz                 Hyx                 Hyy                 Hyz                 Hzx                 Hzy                 Hzz
   1         1        1.2345678901e+00    2.3456789012e-01    3.4567890123e-02    2.3456789012e-01    4.5678901234e+00    5.6789012345e-01    3.4567890123e-02    5.6789012345e-01    6.7890123456e+00
   1         2        -3.4567890123e-01   -4.5678901234e-02   -5.6789012345e-03   -4.5678901234e-02   -6.7890123456e-01   -7.8901234567e-02   -5.6789012345e-03   -7.8901234567e-02   -8.9012345678e-01
   ...

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
       "intervalStep": 10000,
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
