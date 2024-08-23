Steric
------

The Steric potential is used to model steric interactions between particles. It is typically represented as a power law repulsion, which can be useful for preventing overlap between particles or modeling excluded volume effects.

There are two types of Steric potentials implemented:

Steric6
~~~~~~~

.. math::

    U = \epsilon \left(\frac{\sigma}{r}\right)^6

Steric12
~~~~~~~~

.. math::

    U = \epsilon \left(\frac{\sigma}{r}\right)^{12}

where:

* :math:`\epsilon` is the strength of the interaction
* :math:`\sigma` is the characteristic distance
* :math:`r` is the distance between the particles

----

* **type**: ``Bond2``, ``Steric6`` or ``Steric12``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``epsilon``: ``real``: Strength of the interaction :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`

Example (for Steric6):

.. code-block::

   "steric6Bonds":{
     "type":["Bond2","Steric6"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.0, 1.0],
             [1, 2, 0.8, 1.1],
             [2, 3, 1.2, 0.9]]
   }

Steric6Common_epsilon_sigma
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Steric6 bonds variant with common parameters (``epsilon`` and ``sigma``) for all bonds.

----

* **type**: ``Bond2``, ``Steric6Common_epsilon_sigma``
* **parameters**:

  * ``epsilon``: ``real``: Common strength of the interaction for all bonds :math:`[energy]`
  * ``sigma``: ``real``: Common characteristic distance for all bonds :math:`[distance]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle

Example:

.. code-block::

   "steric6BondsCommonEpsilonSigma":{
     "type":["Bond2","Steric6Common_epsilon_sigma"],
     "parameters":{"epsilon":1.0,
                   "sigma":1.0},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [1, 2],
             [2, 3]]
   }

Note: Similar variants (Steric12 and Steric12Common_epsilon_sigma) are available for the Steric12 potential.
