ConstantTorqueOverCenterOfMass
------------------------------

The ConstantTorqueOverCenterOfMass potential applies a constant torque to a set of particles about their center of mass.

The torque :math:`\boldsymbol{\tau}` is applied to the set of particles, causing a rotational acceleration according to:

.. math::

    \mathbf{I} \cdot \boldsymbol{\alpha} = \boldsymbol{\tau}

where :math:`\mathbf{I}` is the moment of inertia tensor and :math:`\boldsymbol{\alpha}` is the angular acceleration.

----

* **type**: ``Set1``, ``ConstantTorqueOverCenterOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the set of particles
  * ``torque``: ``real3``: Constant torque vector :math:`[torque]`

Example:

.. code-block::

   "constantTorque":{
     "type":["Set1","ConstantTorqueOverCenterOfMass"],
     "labels":["idSet_i", "torque"],
     "data":[
       [0, [0.0, 0.0, 1.0]]
     ]
   }

.. note::
   This potential is useful for simulating rotating objects or applying external torques to groups of particles. The torque is distributed among the particles based on their positions relative to the center of mass.

.. warning::
   The energy for this potential is not well-defined, as it represents a non-conservative force.
