SurfaceWCA
----------

The SurfaceWCA potential models the interaction between particles and a planar surface using the Weeks-Chandler-Andersen (WCA) potential. It is a purely repulsive potential, useful for modeling excluded volume effects near surfaces.

There are two types of SurfaceWCA potentials implemented:

SurfaceWCAType1
~~~~~~~~~~~~~~~

.. math::

    U(z) = \begin{cases}
    4\epsilon \left[ \left(\frac{\sigma}{z-z_s}\right)^{12} - \left(\frac{\sigma}{z-z_s}\right)^6 \right] + \epsilon & \text{if } z-z_s < 2^{1/6}\sigma \\
    0 & \text{if } z-z_s \geq 2^{1/6}\sigma
    \end{cases}

SurfaceWCAType2
~~~~~~~~~~~~~~~

.. math::

    U(z) = \begin{cases}
    \epsilon \left[ \left(\frac{\sigma}{z-z_s}\right)^{12} - 2\left(\frac{\sigma}{z-z_s}\right)^6 \right] + \epsilon & \text{if } z-z_s < \sigma \\
    0 & \text{if } z-z_s \geq \sigma
    \end{cases}

where:

* :math:`\epsilon` is the energy scale of the interaction
* :math:`\sigma` is the characteristic distance
* :math:`z` is the distance of the particle from the surface
* :math:`z_s` is the position of the surface

----

* **type**: ``Surface``, ``SurfaceWCAType1`` or ``SurfaceWCAType2``
* **parameters**:

  * ``surfacePosition``: ``real``: Position of the surface along the z-axis :math:`[distance]`

* **data**:

  * ``name``: ``string``: Name of the particle type
  * ``epsilon``: ``real``: Energy scale of the interaction :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`

Example:

.. code-block::

   "surfaceWCA":{
     "type":["Surface","SurfaceWCAType1"],
     "parameters":{
       "surfacePosition": 0.0
     },
     "labels":["name", "epsilon", "sigma"],
     "data":[
       ["A", 1.0, 1.0],
       ["B", 0.5, 0.8]
     ]
   }


.. tip::
   This potential is particularly useful for simulating non-adsorbing surfaces or creating effective walls in simulations.
