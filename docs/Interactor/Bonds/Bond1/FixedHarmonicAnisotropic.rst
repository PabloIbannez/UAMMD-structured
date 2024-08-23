FixedHarmonicAnisotropic
------------------------

The FixedHarmonicAnisotropic potential applies an anisotropic harmonic restraint between a particle and a fixed point in space, allowing for different spring constants in each dimension.

.. math::

    U = \frac{1}{2}K_x(x - x_0)^2 + \frac{1}{2}K_y(y - y_0)^2 + \frac{1}{2}K_z(z - z_0)^2

where:

* :math:`K_x, K_y, K_z` are the spring constants in each dimension
* :math:`x, y, z` are the coordinates of the particle
* :math:`x_0, y_0, z_0` are the coordinates of the fixed point

----

* **type**: ``Bond1``, ``FixedHarmonicAnisotropic``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the particle
  * ``K``: ``real3``: Spring constants in each dimension :math:`[energy]/[distance]^2`
  * ``r0``: ``real3``: Equilibrium distances in each dimension :math:`[distance]`
  * ``position``: ``real3``: Fixed point coordinates :math:`[distance]`

Example:

.. code-block::

   "fixedHarmonicAnisotropicBonds":{
     "type":["Bond1","FixedHarmonicAnisotropic"],
     "parameters":{},
     "labels":["id_i", "K", "r0", "position"],
     "data":[[0, [100.0, 50.0, 75.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]],
             [1, [75.0, 100.0, 50.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]
   }
