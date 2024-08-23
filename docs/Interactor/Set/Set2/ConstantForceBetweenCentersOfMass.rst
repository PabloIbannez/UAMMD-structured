ConstantForceBetweenCentersOfMass
---------------------------------

The ConstantForceBetweenCentersOfMass potential applies a constant force between the centers of mass of two sets of particles.

.. math::

    \mathbf{F} = F\hat{\mathbf{r}}

where:

* :math:`F` is the magnitude of the force
* :math:`\hat{\mathbf{r}}` is the unit vector pointing from the center of mass of set i to the center of mass of set j

The energy of this potential is:

.. math::

    U = -F \cdot r

where :math:`r` is the distance between the centers of mass.

----

* **type**: ``Set2``, ``ConstantForceBetweenCentersOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the first set of particles
  * ``idSet_j``: ``int``: Identifier for the second set of particles
  * ``force``: ``real``: Magnitude of the constant force :math:`[force]`

Example:

.. code-block::

   "constantForce":{
     "type":["Set2","ConstantForceBetweenCentersOfMass"],
     "labels":["idSet_i", "idSet_j", "force"],
     "data":[
       [0, 1, 10.0]
     ]
   }

.. note::
   This potential can be used to model constant attractive or repulsive forces between groups of particles, such as electrostatic interactions between charged objects.
