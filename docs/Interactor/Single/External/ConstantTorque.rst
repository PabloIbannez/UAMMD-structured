ConstantTorque
--------------

The ConstantTorque potential applies a uniform, constant torque to all particles in the system (or a selected group of particles).
This can be used to model external rotational fields or to apply specific torques to selected particles.

----

* **type**: ``External``, ``ConstantTorque``
* **parameters**:

  * ``constantTorque``: ``real3``: The constant torque vector to be applied :math:`[torque]`

Example:

.. code-block::

   "constantTorque":{
     "type":["External","ConstantTorque"],
     "parameters":{
       "constantTorque":[0.0, 0.0, 1.0]
     }
   }

.. note::
   This potential applies the same torque to all particles in the selected group. To apply different torques to different particles, you may need to use multiple ConstantTorque potentials with different particle groups.

.. tip::
   The ConstantTorque potential can be useful for simulating particles in rotating fields, studying the alignment of anisotropic particles, or modeling molecular motors.
