Morse
-----

The Morse potential is an interatomic interaction model that better describes the anharmonicity of real bonds compared to the harmonic potential. It allows bond breaking and is often used for modeling diatomic molecules or metal surfaces.

.. math::

    U = E \left[ 1 - e^{-(r-r0)/D} \right]^2 - E

where:

* :math:`E` is the well depth (depth of the potential energy minimum)
* :math:`D` is a parameter controlling the width of the potential
* :math:`r` is the current distance between the bonded particles
* :math:`r0` is the equilibrium bond distance

----

* **type**: ``Bond2``, ``Morse``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``r0``: ``real``: Equilibrium bond distance :math:`[distance]`
  * ``E``: ``real``: Well depth :math:`[energy]`
  * ``D``: ``real``: Bond width :math:`[distance]`

Example:

.. code-block::

   "morseBonds":{
     "type":["Bond2","Morse"],
     "parameters":{},
     "labels":["id_i", "id_j", "r0", "E", "D"],
     "data":[[0, 1, 1.0, 100.0, 10.0],
             [1, 2, 1.1, 100.0, 10.0],
             [2, 3, 0.9, 100.0, 10.0]]
   }

MorseCommon_D
~~~~~~~~~~~~~

Morse bonds variant with a common :math:`D` parameter for all bonds.

----

* **type**: ``Bond2``, ``MorseCommon_D``
* **parameters**:

  * ``D``: ``real``: Potential width, common for all bonds :math:`[distance]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``r0``: ``real``: Equilibrium bond distance :math:`[distance]`
  * ``E``: ``real``: Well depth :math:`[energy]`

Example:

.. code-block::

   "morseBondsCommonD":{
     "type":["Bond2","MorseCommon_D"],
     "parameters":{"D":10.0},
     "labels":["id_i", "id_j", "r0", "E"],
     "data":[[0, 1, 1.0, 100.0],
             [1, 2, 1.1, 100.0],
             [2, 3, 0.9, 100.0]]
   }

MorseCommon_r0_E_D
~~~~~~~~~~~~~~~~~~

Morse bonds variant with common parameters (``r0``, ``E``, and ``D``) for all bonds.

----

* **type**: ``Bond2``, ``MorseCommon_r0_E_D``
* **parameters**:

  * ``r0``: ``real``: Common equilibrium bond distance for all bonds :math:`[distance]`
  * ``E``: ``real``: Common well depth for all bonds :math:`[energy]`
  * ``D``: ``real``: Common potential width for all bonds :math:`[distance]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle

Example:

.. code-block::

   "morseBondsCommonR0ED":{
     "type":["Bond2","MorseCommon_r0_E_D"],
     "parameters":{"r0":1.0,
                   "E":100.0,
                   "D":10.0},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [1, 2],
             [2, 3]]
   }
