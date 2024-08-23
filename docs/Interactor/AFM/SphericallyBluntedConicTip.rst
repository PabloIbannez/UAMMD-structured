SphericallyBluntedConicTip
--------------------------

The SphericallyBluntedConicTip potential models an AFM tip as a conic structure with a spherical end. This geometry is often used to represent real AFM tips more accurately than simple spherical models.

The potential includes:
1. Tip-sample interaction: A repulsive potential between the tip and sample particles.
2. Tip-chip interaction: A harmonic potential representing the cantilever.

----

* **type**: ``AFM``, ``SphericallyBluntedConicTip``
* **parameters**:
  * ``sigma``: ``real``: Interaction range parameter :math:`[distance]`
  * ``epsilon``: ``real``: Interaction strength parameter :math:`[energy]`
  * ``K``: ``real``: Vertical spring constant of the cantilever :math:`[energy/distance^2]`
  * ``Kxy``: ``real``: Lateral spring constant of the cantilever :math:`[energy/distance^2]`
  * ``tipAngle``: ``real``: Half-angle of the conical part of the tip :math:`[angle]`
  * ``tipVelocity``: ``real``: Velocity of the tip in the z-direction :math:`[distance/time]`
  * ``startChipPosition``: ``real3``: Initial position of the chip (cantilever base) :math:`[distance]`

Example:

.. code-block::

   "afmTip":{
     "type":["AFM","SphericallyBluntedConicTip"],
     "parameters":{
       "sigma": 0.1,
       "epsilon": 100.0,
       "K": 10.0,
       "Kxy": 1.0,
       "tipAngle": 0.261799, // 15 degrees in radians
       "tipVelocity": 0.1,
       "startChipPosition": [0.0, 0.0, 20.0]
     }
   }

.. note::
   This potential is particularly useful for modeling AFM tips with a more realistic geometry, allowing for accurate simulation of both normal and lateral forces during AFM experiments.

.. warning::
   The tip angle should be specified in radians. Ensure that the tip velocity and start position are consistent with your simulation time and length scales.
