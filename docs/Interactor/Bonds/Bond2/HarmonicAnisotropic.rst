HarmonicAnisotropic
-------------------

The HarmonicAnisotropic potential is an extension of the standard Harmonic potential that allows for different spring constants along different spatial directions. This is useful for modeling bonds with directional preferences or in anisotropic environments.

.. math::

    U = \frac{1}{2}(K_x(x-x_0)^2 + K_y(y-y_0)^2 + K_z(z-z_0)^2)

where:

* :math:`K_x, K_y, K_z` are the spring constants in x, y, and z directions
* :math:`x, y, z` are the current component distances between the bonded particles
* :math:`x_0, y_0, z_0` are the equilibrium component distances

----

* **type**: ``Bond2``, ``HarmonicAnisotropic``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``K``   : ``real3``: Spring constants :math:`(K_x, K_y, K_z)` :math:`[energy]/[distance]^2`
  * ``r0``  : ``real3``: Equilibrium distances :math:`(x_0, y_0, z_0)` :math:`[distance]`

Example:

.. code-block::

   "harmonicAnisotropicBonds":{
     "type":["Bond2","HarmonicAnisotropic"],
     "parameters":{},
     "labels":["id_i", "id_j", "K", "r0"],
     "data":[[0, 1, [100.0, 50.0, 75.0], [1.0, 0.5, 0.8]],
             [1, 2, [120.0, 60.0, 90.0], [0.9, 0.6, 0.7]],
             [2, 3, [80.0, 40.0, 60.0], [1.1, 0.4, 0.9]]]
   }

