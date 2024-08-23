ConstantTorqueBetweenCentersOfMass
----------------------------------

The ConstantTorqueBetweenCentersOfMass potential applies a constant torque between two sets of particles, acting on the line connecting their centers of mass.

The torque :math:`\boldsymbol{\tau}` is applied perpendicular to the line connecting the centers of mass, causing a rotational acceleration according to:

.. math::

    \mathbf{I} \cdot \boldsymbol{\alpha} = \boldsymbol{\tau}

where :math:`\mathbf{I}` is the moment of inertia tensor and :math:`\boldsymbol{\alpha}` is the angular acceleration.

----

* **type**: ``Set2``, ``ConstantTorqueBetweenCentersOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the first set of particles
  * ``idSet_j``: ``int``: Identifier for the second set of particles
  * ``torque``: ``real``: Magnitude of the constant torque :math:`[torque]`

Example:

.. code-block::

   "constantTorque":{
     "type":["Set2","ConstantTorqueBetweenCentersOfMass"],
     "labels":["idSet_i", "idSet_j", "torque"],
     "data":[
       [0, 1, 1.0]
     ]
   }

.. note::
   This potential is useful for simulating rotational forces between two groups of particles, such as those found in molecular motors or other complex mechanical systems.

.. warning::
   The energy for this potential is not well-defined, as it represents a non-conservative force. The torque is distributed among the particles in each set based on their positions relative to their respective centers of mass.
