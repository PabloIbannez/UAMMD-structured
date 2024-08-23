ConstantForce
-------------

The ConstantForce potential applies a constant force to a particle. This can be used to model external fields or to apply specific forces to selected particles.

.. math::

    U = -\mathbf{F}\cdot \mathbf{r}

----

* **type**: ``Bond1``, ``ConstantForce``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the particle
  * ``force``: ``real3``: Constant force vector :math:`[energy]/[distance]`

Example:

.. code-block::

   "constantForceBonds":{
     "type":["Bond1","ConstantForce"],
     "parameters":{},
     "labels":["id_i", "force"],
     "data":[[0, [1.0, 0.0, 0.0]],
             [1, [0.0, -1.0, 0.5]]]
   }
