SphericalTip
------------

The SphericalTip potential models an AFM tip as a sphere. This is a simpler model compared to the SphericallyBluntedConicTip, but it can be sufficient for many applications and is computationally less expensive.

The potential includes:
1. Tip-sample interaction: A Lennard-Jones-like potential between the tip and sample particles.
2. Tip-chip interaction: A harmonic potential representing the cantilever.

This potential also includes features for simulating approach and retraction curves in AFM indentation experiments.

----

* **type**: ``AFM``, ``SphericalTip``
* **parameters**:
  * ``epsilon``: ``real``: Interaction strength parameter :math:`[energy]`
  * ``sigma``: ``real``: Interaction range parameter :math:`[distance]`
  * ``K``: ``real``: Vertical spring constant of the cantilever :math:`[energy/distance^2]`
  * ``Kxy``: ``real``: Lateral spring constant of the cantilever :math:`[energy/distance^2]`
  * ``tipVelocity``: ``real``: Velocity of the tip in the z-direction :math:`[distance/time]`
  * ``startChipPosition``: ``real3``: Initial position of the chip (cantilever base) :math:`[distance]`
  * ``indentationStartStep``: ``ullint``: Simulation step to start the indentation
  * ``indentationBackwardStep``: ``ullint``: Simulation step to start the retraction (optional)

Example:

.. code-block::

   "afmTip":{
     "type":["AFM","SphericalTip"],
     "parameters":{
       "epsilon": 100.0,
       "sigma": 0.1,
       "K": 10.0,
       "Kxy": 1.0,
       "tipVelocity": 0.1,
       "startChipPosition": [0.0, 0.0, 20.0],
       "indentationStartStep": 1000,
       "indentationBackwardStep": 10000
     }
   }

.. note::
   This potential allows for simulating complete AFM indentation cycles, including approach, contact, and retraction phases. The indentation process starts at ``indentationStartStep`` and reverses at ``indentationBackwardStep`` if specified.

.. tip::
   If ``indentationBackwardStep`` is not specified or set to 0, the tip will continue to indent without retracting.

.. warning::
   Ensure that ``indentationBackwardStep`` is larger than ``indentationStartStep`` if both are specified.
