LennardJonesSoftCore
--------------------

The LennardJonesSoftCore potential is a modified version of the Lennard-Jones potential that introduces a "soft core" to avoid the singularity at r = 0. This is particularly useful in free energy calculations and simulations involving the creation or annihilation of particles.

There are two types of LennardJonesSoftCore potentials implemented:

LennardJonesSoftCoreType1
~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U = 4\epsilon \lambda^n \left[ \frac{1}{(\alpha(1-\lambda)^2 + (r/\sigma)^6)^2} - \frac{1}{\alpha(1-\lambda)^2 + (r/\sigma)^6} \right]

LennardJonesSoftCoreType2
~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U = \epsilon \lambda^n \left[ \frac{1}{(\alpha(1-\lambda)^2 + (r/\sigma)^6)^2} - 2\frac{1}{\alpha(1-\lambda)^2 + (r/\sigma)^6} \right]

where:

* :math:`\epsilon` is the depth of the potential well
* :math:`\sigma` is the distance at which the potential is zero
* :math:`r` is the distance between the particles
* :math:`\lambda` is the coupling parameter (0 ≤ λ ≤ 1)
* :math:`\alpha` is the soft-core parameter
* :math:`n` is the λ-exponent

----

* **type**: ``Bond2``, ``LennardJonesSoftCoreType1`` or ``LennardJonesSoftCoreType2``
* **parameters**:

  * ``alpha``: ``real``: Soft-core parameter
  * ``n``: ``int``: λ-exponent (default: 2)

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``epsilon``: ``real``: Depth of the potential well :math:`[energy]`
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example (for LennardJonesSoftCoreType1):

.. code-block::

   "ljSoftCoreBonds":{
     "type":["Bond2","LennardJonesSoftCoreType1"],
     "parameters":{"alpha":0.5,
                   "n":2},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.0, 1.0],
             [1, 2, 0.8, 1.1]]
   }

LennardJonesSoftCoreType1Common_epsilon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LennardJonesSoftCoreType1 bonds variant with a common well depth (``epsilon``) for all bonds.

----

* **type**: ``Bond2``, ``LennardJonesSoftCoreType1Common_epsilon``
* **parameters**:

  * ``alpha``: ``real``: Soft-core parameter
  * ``n``: ``int``: λ-exponent (default: 2)
  * ``epsilon``: ``real``: Common depth of the potential well for all bonds :math:`[energy]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example:

.. code-block::

   "ljSoftCoreBondsCommonEpsilon":{
     "type":["Bond2","LennardJonesSoftCoreType1Common_epsilon"],
     "parameters":{"alpha":0.5,
                   "n":2,
                   "epsilon":1.0},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 1.0],
             [1, 2, 1.1]]
   }

Note: Similar variants (LennardJonesSoftCoreType2Common_epsilon) are available for LennardJonesSoftCoreType2.
