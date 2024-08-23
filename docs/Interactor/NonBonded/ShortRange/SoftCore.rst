SoftCore
--------

The SoftCore potential is a modified version of the Lennard-Jones potential that removes the singularity at r = 0. It's particularly useful in free energy calculations and simulations involving the creation or annihilation of particles.

There are two types of SoftCore potentials implemented:

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
* :math:`r` is the distance between particles
* :math:`\lambda` is the coupling parameter (0 ≤ :math:`\lambda` ≤ 1)
* :math:`\alpha` is the soft-core parameter
* :math:`n` is the :math:`\lambda`-exponent

----

* **type**: ``NonBonded``, ``LennardJonesSoftCoreType1`` or ``LennardJonesSoftCoreType2``
* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma
  * ``alpha``: ``real``: Soft-core parameter
  * ``n``: ``int``: :math:`\lambda`-exponent (default: 2)

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Depth of the potential well :math:`[energy]`
  * ``sigma``: ``real``: Distance at which the potential is zero :math:`[distance]`

Example:

.. code-block::

   "softCore":{
     "type":["NonBonded","LennardJonesSoftCoreType1"],
     "parameters":{
       "cutOffFactor":2.5,
       "alpha":0.5,
       "n":2,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 1.2, 0.9],
       ["B", "B", 0.8, 1.1]
     ]
   }

.. note::
   The SoftCore potential is particularly useful in alchemical free energy calculations, where particles are gradually introduced or removed from the system. The :math:`\lambda` parameter is typically controlled globally for the entire simulation and is not specified per interaction pair.
