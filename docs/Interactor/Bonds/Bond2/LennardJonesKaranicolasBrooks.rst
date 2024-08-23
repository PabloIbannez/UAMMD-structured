LennardJonesKaranicolasBrooks
------------------------------

The LennardJonesKaranicolasBrooks potential is a modified version of the Lennard-Jones potential, introduced by Karanicolas and Brooks to model protein folding [Karanicolas2002]_.
It includes a well depth scaling factor that depends on the types of interacting residues.

.. math::

    U = \epsilon_{ij} \left[ 13\left(\frac{\sigma_{ij}}{r}\right)^{12} - 18\left(\frac{\sigma_{ij}}{r}\right)^{10} + 4\left(\frac{\sigma_{ij}}{r}\right)^6 \right]

where:

* :math:`\epsilon_{ij}` is the interaction energy between the particles
* :math:`\sigma_{ij}` is the distance at which the potential is zero
* :math:`r` is the distance between the particles

----

* **type**: ``Bond2``, ``LennardJonesKaranicolasBrooks``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``epsilon``: ``real``: Well depth :math:`[energy]`
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example:

.. code-block::

   "ljKBBonds":{
     "type":["Bond2","LennardJonesKaranicolasBrooks"],
     "parameters":{},
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.0, 1.0],
             [1, 2, 0.8, 1.1]]
   }

LennardJonesKaranicolasBrooksCommon_epsilon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LennardJonesKaranicolasBrooks bonds variant with a common well depth (``epsilon``) for all bonds.

----

* **type**: ``Bond2``, ``LennardJonesKaranicolasBrooksCommon_epsilon``
* **parameters**:

  * ``epsilon``: ``real``: Common well depth for all bonds :math:`[energy]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example:

.. code-block::

   "ljKBBondsCommonEpsilon":{
     "type":["Bond2","LennardJonesKaranicolasBrooksCommon_epsilon"],
     "parameters":{"epsilon":1.0},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 1.0],
             [1, 2, 1.1]]
   }

.. [Karanicolas2002] Karanicolas, John, and Charles L. Brooks III. "The origins of asymmetry in the folding transition states of protein L and protein G." Protein Science 11, no. 10 (2002): 2351-2361.
