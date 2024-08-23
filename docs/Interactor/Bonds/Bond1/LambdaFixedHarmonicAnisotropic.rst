LambdaFixedHarmonicAnisotropic
------------------------------

The LambdaFixedHarmonicAnisotropic potential is a variant of the FixedHarmonicAnisotropic potential that includes a coupling parameter :math:`\lambda`, useful for free energy calculations or gradual introduction of restraints.

.. math::

    U = \lambda^n \left( \frac{1}{2}K_x(x - x_0)^2 + \frac{1}{2}K_y(y - y_0)^2 + \frac{1}{2}K_z(z - z_0)^2 \right)

where:

* :math:`\lambda` is the coupling parameter (0 ≤ :math:`\lambda` ≤ 1)
* :math:`n` is the :math:`\lambda`-exponent
* :math:`K_x, K_y, K_z` are the spring constants in each dimension
* :math:`x, y, z` are the coordinates of the particle
* :math:`x_0, y_0, z_0` are the coordinates of the fixed point

----

* **type**: ``Bond1``, ``LambdaFixedHarmonicAnisotropic``
* **parameters**:

  * ``n``: ``int``: :math:`\lambda`-exponent (default: 2)

* **data**:

  * ``id_i``: ``int``: Id of the particle
  * ``K``: ``real3``: Spring constants in each dimension :math:`[energy]/[distance]^2`
  * ``r0``: ``real3``: Equilibrium distances in each dimension :math:`[distance]`
  * ``position``: ``real3``: Fixed point coordinates :math:`[distance]`

Example:

.. code-block::

   "lambdaFixedHarmonicAnisotropicBonds":{
     "type":["Bond1","LambdaFixedHarmonicAnisotropic"],
     "parameters":{"n":2},
     "labels":["id_i", "K", "r0", "position"],
     "data":[[0, [100.0, 50.0, 75.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]],
             [1, [75.0, 100.0, 50.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]
   }

.. note::
    The :math:`\lambda` parameter is typically controlled globally for the entire simulation and is not specified per bond.
