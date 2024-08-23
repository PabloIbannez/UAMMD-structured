FixedHarmonicAnisotropicCenterOfMass
------------------------------------

The FixedHarmonicAnisotropicCenterOfMass potential applies an anisotropic harmonic restraint to the center of mass of a set of particles, tethering it to a fixed point in space.

.. math::

    U = \frac{1}{2}K_x(x-x_0)^2 + \frac{1}{2}K_y(y-y_0)^2 + \frac{1}{2}K_z(z-z_0)^2

where:

* :math:`K_x, K_y, K_z` are the spring constants in each direction
* :math:`x, y, z` are the coordinates of the center of mass
* :math:`x_0, y_0, z_0` are the coordinates of the fixed point

----

* **type**: ``Set1``, ``FixedHarmonicAnisotropicCenterOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the set of particles
  * ``K``: ``real3``: Spring constants in x, y, and z directions :math:`[energy/distance^2]`
  * ``r0``: ``real3``: Equilibrium position :math:`[distance]`
  * ``position``: ``real3``: Fixed point coordinates :math:`[distance]`

Example:

.. code-block::

   "fixedHarmonicAnisotropic":{
     "type":["Set1","FixedHarmonicAnisotropicCenterOfMass"],
     "labels":["idSet_i", "K", "r0", "position"],
     "data":[
       [0, [100.0, 100.0, 100.0], [0.0, 0.0, 0.0], [5.0, 5.0, 5.0]]
     ]
   }

.. note::
   This potential is useful for creating flexible tethers or anchors for groups of particles, allowing anisotropic behavior in different directions.
