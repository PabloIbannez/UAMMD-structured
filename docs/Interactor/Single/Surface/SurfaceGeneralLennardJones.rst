SurfaceGeneralLennardJones
--------------------------

The SurfaceGeneralLennardJones potential is a flexible surface interaction that combines properties of both the Lennard-Jones and WCA potentials. It allows for both attractive and purely repulsive interactions based on the sign of the epsilon parameter.

There are two types of SurfaceGeneralLennardJones potentials implemented:

SurfaceGeneralLennardJonesType1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U(z) = \begin{cases}
    WCA_{Type1}(z-z_s; \epsilon, \sigma) & \text{if } \epsilon > 0 \\
    LJ_{Type1}(z-z_s; |\epsilon|, \sigma) & \text{if } \epsilon < 0
    \end{cases}

SurfaceGeneralLennardJonesType2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U(z) = \begin{cases}
    WCA_{Type2}(z-z_s; \epsilon, \sigma) & \text{if } \epsilon > 0 \\
    LJ_{Type2}(z-z_s; |\epsilon|, \sigma) & \text{if } \epsilon < 0
    \end{cases}

where :math:`WCA` and :math:`LJ` refer to the respective WCA and Lennard-Jones potentials defined in the SurfaceWCA and SurfaceLennardJones sections.

----

* **type**: ``Surface``, ``SurfaceGeneralLennardJonesType1`` or ``SurfaceGeneralLennardJonesType2``
* **parameters**:

  * ``surfacePosition``: ``real``: Position of the surface along the z-axis :math:`[distance]`

* **data**:

  * ``name``: ``string``: Name of the particle type
  * ``epsilon``: ``real``: Energy scale of the interaction (sign determines attractive/repulsive nature) :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`

Example:

.. code-block::

   "surfaceGeneralLJ":{
     "type":["Surface","SurfaceGeneralLennardJonesType1"],
     "parameters":{
       "surfacePosition": 0.0
     },
     "labels":["name", "epsilon", "sigma"],
     "data":[
       ["A", 1.0, 1.0],
       ["B", -0.5, 0.8]
     ]
   }

.. note::
   The SurfaceGeneralLennardJones potential behaves like WCA for positive epsilon values and like Lennard-Jones for negative epsilon values.

.. tip::
   This potential is versatile and can be used to model a wide range of surface interactions, from purely repulsive to strongly attractive, within a single simulation.
