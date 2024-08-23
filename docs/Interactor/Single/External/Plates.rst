Plates
------

The Plates potential confines particles between two parallel plates, applying repulsive forces when particles approach the plates.

The initial separation between the plates is set by the ``platesSeparation`` parameter. 
The separation can be changed during the simulation by setting the ``compressionVelocity`` parameter to a non-zero value. 
If ``compressionVelocity`` is positive, the plates will move towards each other at the specified rate, if negative, they will move apart.
The minimum and maximum allowed plate separations can be set using the ``minPlatesSeparation`` and ``maxPlatesSeparation`` parameters, respectively.

The potential uses a steep repulsive interaction (r^12) for each plate:

.. math::

   U(z) = \epsilon \left[\left(\frac{\sigma}{z - z_{\text{bottom}}}\right)^{12} + \left(\frac{\sigma}{z_{\text{top}} - z}\right)^{12}\right]

where :math:`z` is the z-coordinate of the particle, :math:`z_{\text{bottom}}` and :math:`z_{\text{top}}` are the positions of the bottom and top plates, :math:`\epsilon` is the energy scale, and :math:`\sigma` is the length scale of the interaction.


----

* **type**: ``External``, ``Plates``
* **parameters**:

  * ``platesSeparation``: ``real``: Initial separation between the plates :math:`[distance]`
  * ``platesEpsilon``: ``real``: Energy scale of the repulsive interaction :math:`[energy]` (default: 1.0)
  * ``platesSigma``: ``real``: Length scale of the repulsive interaction :math:`[distance]` (default: 1.0)
  * ``compressionVelocity``: ``real``: Rate of change of the plate separation :math:`[distance/time]` (default: 0.0)
  * ``minPlatesSeparation``: ``real``: Minimum allowed plate separation :math:`[distance]` (default: 0.0)
  * ``maxPlatesSeparation``: ``real``: Maximum allowed plate separation :math:`[distance]` (default: infinity)

Example:

.. code-block::

   "plates":{
     "type":["External","Plates"],
     "parameters":{
       "platesSeparation": 20.0,
       "platesEpsilon": 1.0,
       "platesSigma": 1.0,
       "compressionVelocity": 0.01,
       "minPlatesSeparation": 10.0,
       "maxPlatesSeparation": 30.0
     }
   }

.. warning::
   Ensure that the initial configuration of particles is between the plates to avoid extremely large forces at the start of the simulation.

.. tip::
   This potential is useful for simulating confined systems or studying the behavior of materials under compression.
