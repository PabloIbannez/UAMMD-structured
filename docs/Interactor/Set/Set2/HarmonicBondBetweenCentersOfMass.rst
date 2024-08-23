HarmonicBondBetweenCentersOfMass
--------------------------------

The HarmonicBondBetweenCentersOfMass potential applies a harmonic bond between the centers of mass of two sets of particles.

.. math::

    U = \frac{1}{2}K(r-r_0)^2

where:

* :math:`K` is the spring constant
* :math:`r` is the distance between the centers of mass of the two sets
* :math:`r_0` is the equilibrium distance

----

* **type**: ``Set2``, ``HarmonicBondBetweenCentersOfMass``
* **parameters**: None
* **data**:
  * ``idSet_i``: ``int``: Identifier for the first set of particles
  * ``idSet_j``: ``int``: Identifier for the second set of particles
  * ``K``: ``real``: Spring constant :math:`[energy/distance^2]`
  * ``r0``: ``real``: Equilibrium distance :math:`[distance]`

Example:

.. code-block::

   "harmonicBond":{
     "type":["Set2","HarmonicBondBetweenCentersOfMass"],
     "labels":["idSet_i", "idSet_j", "K", "r0"],
     "data":[
       [0, 1, 100.0, 5.0]
     ]
   }

.. note::
   This potential is useful for creating flexible connections between groups of particles or for maintaining a preferred distance between two structures in a simulation.
