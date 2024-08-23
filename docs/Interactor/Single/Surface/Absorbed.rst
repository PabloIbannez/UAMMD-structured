Absorbed
--------

The Absorbed potential models particles that have been absorbed onto a surface, applying a harmonic restraint to keep them near their absorption point. This is useful for simulating irreversible adsorption or creating stable surface-bound particles.

.. math::

    U(\mathbf{r}) = \begin{cases}
    \frac{1}{2}K((x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2) & \text{if particle is absorbed} \\
    0 & \text{otherwise}
    \end{cases}

where:

* :math:`K` is the spring constant of the harmonic restraint
* :math:`\mathbf{r} = (x, y, z)` is the position of the particle
* :math:`\mathbf{r}_0 = (x_0, y_0, z_0)` is the absorption point

The Absorbed potential is applied only to particles that were below the `heightThreshold` when the potential was first computed. 

Note the strategy, playing with the common parameter ``startStep`` we can control the time when the particles are absorbed.
For example we can wait for a structure to relax before it is absorbed. 

----

* **type**: ``Surface``, ``Absorbed``
* **parameters**:

  * ``K``: ``real``: Spring constant of the harmonic restraint :math:`[energy/distance^2]`
  * ``heightThreshold``: ``real``: Height below which particles are considered absorbed :math:`[distance]`

Example:

.. code-block::

   "absorbed":{
     "type":["Surface","Absorbed"],
     "parameters":{
       "K": 100.0,
       "heightThreshold": 1.0
     }
   }

.. note::

   The Absorbed potential is applied only to particles that were below the `heightThreshold` when the potential was first computed. The absorption points are determined at this time and remain fixed throughout the simulation.

.. warning::

   This potential determines the absorbed particles only once, when it is first computed. Ensure that your initial configuration is consistent with your intended absorption state.

.. tip::

   The Absorbed potential can be used to model particles that are irreversibly bound to a surface while still allowing for some movement around their binding point.
