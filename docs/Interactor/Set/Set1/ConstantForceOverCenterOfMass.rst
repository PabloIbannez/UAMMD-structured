ConstantForceOverCenterOfMass
-----------------------------

The ConstantForceOverCenterOfMass potential applies a constant force to the center of mass of a set of particles.

.. math::

    \mathbf{F} = \mathbf{F}_\text{constant}

The energy of this potential is zero, as it represents a non-conservative force.

----

* **type**: ``Set1``, ``ConstantForceOverCenterOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the set of particles
  * ``force``: ``real3``: Constant force vector :math:`[force]`

Example:

.. code-block::

   "constantForce":{
     "type":["Set1","ConstantForceOverCenterOfMass"],
     "labels":["idSet_i", "force"],
     "data":[
       [0, [0.0, 0.0, -9.8]]
     ]
   }

.. note::
   This potential can be used to model external forces acting on entire groups of particles, such as gravity or electric fields.
