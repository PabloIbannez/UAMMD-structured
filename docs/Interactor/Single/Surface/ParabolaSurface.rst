ParabolaSurface
---------------

The ParabolaSurface potential models the interaction between particles and a planar surface using a parabolic potential. It acts only in the z-direction and is one-sided, applying a force only when particles are below the surface.

.. math::

    U(z) = \begin{cases}
    \frac{1}{2}\epsilon(z - z_s)^2, & \text{if } z \leq z_s \\
    0, & \text{if } z > z_s
    \end{cases}

where:

* :math:`\epsilon` is the strength of the potential
* :math:`z_s` is the position of the surface
* :math:`z` is the z-coordinate of the particle

----

* **type**: ``Surface``, ``ParabolaSurface``
* **parameters**:

  * ``surfacePosition``: ``real``: Position of the surface (z_s) :math:`[distance]`
* **data**:

  * ``name``: ``string``: Name of the particle type
  * ``epsilon``: ``real``: Strength of the potential :math:`[energy/distance^2]`
  * ``sigma``: ``real``: Not used, can be set to any value

Example:

.. code-block::

   "parabolaSurface":{
     "type":["Surface","ParabolaSurface"],
     "parameters":{
       "surfacePosition": 0.0
     },
     "labels":["name", "epsilon", "sigma"],
     "data":[
       ["A", 10.0, 1.0],
       ["B", 5.0, 1.0]
     ]
   }

.. note::

   The ParabolaSurface potential only applies forces in the z-direction and only when particles are below the surface (z ≤ z_s). It does not affect particles above the surface.

.. tip::

   This potential can be used to create a soft, one-sided wall or to model interactions with a planar surface where the repulsion increases quadratically with distance.
