SurfaceLennardJones
-------------------

The SurfaceLennardJones potential models the interaction between particles and a planar surface using the Lennard-Jones potential. It is useful for simulating surface adsorption, particle-wall interactions, and other surface phenomena.

There are two types of SurfaceLennardJones potentials implemented:

SurfaceLennardJonesType1
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U(z) = 4\epsilon \left[ \left(\frac{\sigma}{z-z_s}\right)^{12} - \left(\frac{\sigma}{z-z_s}\right)^6 \right]

SurfaceLennardJonesType2
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U(z) = \epsilon \left[ \left(\frac{\sigma}{z-z_s}\right)^{12} - 2\left(\frac{\sigma}{z-z_s}\right)^6 \right]

where:

* :math:`\epsilon` is the depth of the potential well
* :math:`\sigma` is the distance at which the potential is zero
* :math:`z` is the distance of the particle from the surface
* :math:`z_s` is the position of the surface

----

* **type**: ``Surface``, ``SurfaceLennardJonesType1`` or ``SurfaceLennardJonesType2``
* **parameters**:

  * ``surfacePosition``: ``real``: Position of the surface along the z-axis :math:`[distance]`

* **data**:

  * ``name``: ``string``: Name of the particle type
  * ``epsilon``: ``real``: Depth of the potential well :math:`[energy]`
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example:

.. code-block::

   "surfaceLJ":{
     "type":["Surface","SurfaceLennardJonesType1"],
     "parameters":{
       "surfacePosition": 0.0
     },
     "labels":["name", "epsilon", "sigma"],
     "data":[
       ["A", 1.0, 1.0],
       ["B", 0.5, 0.8]
     ]
   }

.. note::
   The force is only applied in the z-direction, perpendicular to the surface.

