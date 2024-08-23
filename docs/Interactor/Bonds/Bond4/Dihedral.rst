Dihedral
--------

The Dihedral potential is a common form used to model the energetics of dihedral angles in molecular simulations. It is often expressed as a cosine series.

.. math::

    U = K[1 + \cos(n\phi - \phi_0)]

where:

* :math:`K` is the force constant
* :math:`n` is the multiplicity (number of minima as the bond is rotated through 360°)
* :math:`\phi` is the dihedral angle
* :math:`\phi_0` is the phase shift

----

* **type**: ``Bond4``, ``Dihedral``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the second particle
  * ``id_k``: ``int``: Id of the third particle
  * ``id_l``: ``int``: Id of the fourth particle
  * ``n``: ``int``: Multiplicity
  * ``K``: ``real``: Force constant :math:`[energy]`
  * ``phi0``: ``real``: Phase shift :math:`[angle]`

Example:

.. code-block::

   "dihedralBonds":{
     "type":["Bond4","Dihedral"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k", "id_l", "n", "K", "phi0"],
     "data":[[0, 1, 2, 3, 3, 1.0, 0.0],
             [1, 2, 3, 4, 2, 0.5, 3.14]]
   }

DihedralCommon_n_K_phi0
~~~~~~~~~~~~~~~~~~~~~~~

Dihedral bonds variant with common parameters (``n``, ``K``, and ``phi0``) for all bonds.

----

* **type**: ``Bond4``, ``DihedralCommon_n_K_phi0``
* **parameters**:

  * ``n``: ``int``: Common multiplicity for all bonds
  * ``K``: ``real``: Common force constant for all bonds :math:`[energy]`
  * ``phi0``: ``real``: Common phase shift for all bonds :math:`[angle]`

* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the second particle
  * ``id_k``: ``int``: Id of the third particle
  * ``id_l``: ``int``: Id of the fourth particle

Example:

.. code-block::

   "dihedralBondsCommon":{
     "type":["Bond4","DihedralCommon_n_K_phi0"],
     "parameters":{"n":3,
                   "K":1.0,
                   "phi0":0.0},
     "labels":["id_i", "id_j", "id_k", "id_l"],
     "data":[[0, 1, 2, 3],
             [1, 2, 3, 4]]
   }

Dihedral4
~~~~~~~~~

A variant of the Dihedral potential that includes up to four terms in the cosine series.

.. math::

    U = \sum_{n=1}^4 K_n[1 + \cos(n\phi - \phi_{0,n})]

----

* **type**: ``Bond4``, ``Dihedral4``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the second particle
  * ``id_k``: ``int``: Id of the third particle
  * ``id_l``: ``int``: Id of the fourth particle
  * ``K``: ``real4``: Force constants for n=1,2,3,4 :math:`[energy]`
  * ``phi0``: ``real4``: Phase shifts for n=1,2,3,4 :math:`[angle]`

Example:

.. code-block::

   "dihedral4Bonds":{
     "type":["Bond4","Dihedral4"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k", "id_l", "K", "phi0"],
     "data":[[0, 1, 2, 3, [1.0, 0.5, 0.25, 0.1], [0.0, 3.14, 1.57, 0.0]],
             [1, 2, 3, 4, [0.8, 0.4, 0.2, 0.05], [1.57, 0.0, 3.14, 1.57]]]
   }
