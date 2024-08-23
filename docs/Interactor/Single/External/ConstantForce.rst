ConstantForce
-------------

The ConstantForce potential applies a uniform, constant force to all particles in the system (or in a selected group).
This can be used to model external fields or to apply specific forces to selected particles.

----

* **type**: ``External``, ``ConstantForce``
* **parameters**:

  * ``constantForce``: ``real3``: The constant force vector to be applied :math:`[force]`

Example:

.. code-block::

   "constantForce":{
     "type":["External","ConstantForce"],
     "parameters":{
       "constantForce":[0.0, 0.0, -9.8]
     }
   }

.. note::
   The ConstantForce potential applies the same force to all particles in the selected group. To apply different forces to different particles, you may need to use multiple ConstantForce potentials with different particle groups.

.. tip::
   This potential can be useful for simulating gravity, electric fields, or other uniform external forces.
