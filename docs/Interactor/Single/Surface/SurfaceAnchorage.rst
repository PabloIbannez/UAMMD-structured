SurfaceAnchorage
----------------

The SurfaceAnchorage potential models the interaction between particles and a planar surface using a harmonic well within a certain range. It is useful for simulating tethered or anchored particles on a surface.

.. math::

    U(z) = \begin{cases}
    \frac{1}{2}\epsilon\left(\frac{z - z_s}{\sigma}\right)^2 - \frac{1}{2}\epsilon, & \text{if } |z - z_s| \leq \sigma \\
    0, & \text{otherwise}
    \end{cases}

where:

* :math:`\epsilon` is the well depth
* :math:`\sigma` is the width of the well
* :math:`z` is the z-coordinate of the particle
* :math:`z_s` is the position of the surface

----

* **type**: ``Surface``, ``SurfaceAnchorage``
* **parameters**:

  * ``surfacePosition``: ``real``: Position of the surface along the z-axis :math:`[distance]`
* **data**:

  * ``name``: ``string``: Name of the particle type
  * ``epsilon``: ``real``: Well depth :math:`\epsilon` :math:`[energy]`
  * ``sigma``: ``real``: Width of the well :math:`\sigma` :math:`[distance]`

Example:

.. code-block::

   "surfaceAnchorage":{
     "type":["Surface","SurfaceAnchorage"],
     "parameters":{
       "surfacePosition": 0.0
     },
     "labels":["name", "epsilon", "sigma"],
     "data":[
       ["A", 10.0, 1.0],
       ["B", 5.0, 1.5]
     ]
   }

.. note::

   The force is only applied in the z-direction, perpendicular to the surface, and only within the range :math:`|z - z_s| \leq \sigma`. Outside this range, both the energy and force are zero.

.. tip::

   SurfaceAnchorage can be used to model particles tethered to a surface with some flexibility or to create a soft wall boundary condition with a finite range of interaction.
