SphericalShell
--------------

The SphericalShell potential confines particles within a spherical shell, 
applying a repulsive force when particles approach the shell boundary.

The shell is initialized with a center and radius, and particles are confined within the shell by a steep repulsive interaction.
If a non-zero ``radiusVelocity`` is specified, the shell radius will change over time. It will increase if the radiusVelocity is positive, and decrease if it is negative.
It will stop changing when the shell radius reaches the ``minShellRadius`` or ``maxShellRadius``.

The potential uses a steep repulsive interaction (r^12) to confine particles:

.. math::

   U(r) = \epsilon \left(\frac{\sigma}{r - R}\right)^{12}

where :math:`r` is the distance from the particle to the shell center, :math:`R` is the shell radius, :math:`\epsilon` is the energy scale, and :math:`\sigma` is the length scale of the interaction.


----

* **type**: ``External``, ``SphericalShell``
* **parameters**:

  * ``shellCenter``: ``real3``: Center coordinates of the spherical shell
  * ``shellRadius``: ``real``: Radius of the spherical shell :math:`[distance]`
  * ``shellEpsilon``: ``real``: Energy scale of the repulsive interaction :math:`[energy]` (default: 1.0)
  * ``shellSigma``: ``real``: Length scale of the repulsive interaction :math:`[distance]` (default: 1.0)
  * ``minShellRadius``: ``real``: Minimum allowed shell radius :math:`[distance]` (default: 0.0)
  * ``maxShellRadius``: ``real``: Maximum allowed shell radius :math:`[distance]` (default: infinity)
  * ``radiusVelocity``: ``real``: Rate of change of the shell radius :math:`[distance/time]` (default: 0.0)

Example:

.. code-block::

   "sphericalShell":{
     "type":["External","SphericalShell"],
     "parameters":{
       "shellCenter": [0.0, 0.0, 0.0],
       "shellRadius": 10.0,
       "shellEpsilon": 1.0,
       "shellSigma": 1.0,
       "minShellRadius": 5.0,
       "maxShellRadius": 15.0,
       "radiusVelocity": 0.1
     }
   }

.. warning::
   Ensure that the initial configuration of particles is within the spherical shell to avoid extremely large forces at the start of the simulation.
