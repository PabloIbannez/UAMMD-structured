Harmonic
--------

The Harmonic potential is a simple and widely used model for bonded interactions in molecular simulations. It represents a linear spring-like force between two particles, providing a quadratic energy penalty for deviations from the equilibrium bond length.

.. math::

    U = \frac{1}{2}K(r-r_0)^2

where:

- :math:`K` is the spring constant
- :math:`r` is the current distance between the bonded particles
- :math:`r_0` is the equilibrium distance

The Harmonic potential is suitable for modeling covalent bonds in many molecular systems, especially when bond stretching is expected to be small. It's computationally efficient and provides a good approximation for small deviations from the equilibrium length.

----

* **type**: ``Bond2``, ``Harmonic``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``K``   : ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``r0``  : ``real``: Equilibrium distance :math:`[distance]`

Example:

.. code-block::

   "harmonicBonds":{
     "type":["Bond2","Harmonic"],
     "parameters":{},
     "labels":["id_i", "id_j", "K", "r0"],
     "data":[[0, 1, 100.0, 1.0],
             [1, 2, 100.0, 1.0],
             [2, 3, 100.0, 1.0]]
   }

HarmonicCommon_K
~~~~~~~~~~~~~~~~

Harmonic bonds variant with a common spring constant (``K``) for all bonds.

----

* **type**: ``Bond2``, ``HarmonicCommon_K``
* **parameters**:

  * ``K``: ``real``: Common spring constant for all bonds :math:`[energy]/[distance]^2`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``r0``  : ``real``: Equilibrium distance :math:`[distance]`

Example:

.. code-block::

   "harmonicBondsCommonK":{
     "type":["Bond2","HarmonicCommon_K"],
     "parameters":{"K":100.0},
     "labels":["id_i", "id_j", "r0"],
     "data":[[0, 1, 1.0],
             [1, 2, 1.1],
             [2, 3, 0.9]]
   }

HarmonicCommon_K_r0
~~~~~~~~~~~~~~~~~~~

Harmonic bonds variant with common parameters (``K`` and ``r0``) for all bonds.

----

* **type**: ``Bond2``, ``HarmonicCommon_K_r0``
* **parameters**:

  * ``K`` : ``real``: Common spring constant for all bonds :math:`[energy]/[distance]^2`
  * ``r0``: ``real``: Common equilibrium distance for all bonds :math:`[distance]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle

Example:

.. code-block::

   "harmonicBondsCommonKR0":{
     "type":["Bond2","HarmonicCommon_K_r0"],
     "parameters":{"K":100.0,
                   "r0":1.0},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [1, 2],
             [2, 3]]
   }
