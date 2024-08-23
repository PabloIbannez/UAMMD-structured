GeneralLennardJones
-------------------

The GeneralLennardJones potential is an extension of the standard Lennard-Jones potential that allows for both attractive and repulsive interactions based on the sign of the :math:`\epsilon` parameter.
It combines the properties of the Lennard-Jones and Weeks-Chandler-Andersen (WCA) potentials.

There are three types of GeneralLennardJones potentials implemented:

GeneralLennardJonesType1
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U = \begin{cases}
     WCA_{\text{Type1}}(r;\epsilon,\sigma) & \text{if } \epsilon > 0 \\
     LJ_{\text{Type1}}(r;|\epsilon|,\sigma) & \text{if } \epsilon < 0
    \end{cases}

GeneralLennardJonesType2
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U = \begin{cases}
     WCA_{\text{Type2}}(r;\epsilon,\sigma) & \text{if } \epsilon > 0 \\
     LJ_{\text{Type2}}(r;|\epsilon|,\sigma) & \text{if } \epsilon < 0
    \end{cases}

GeneralLennardJonesType3
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    U = \begin{cases}
     WCA_{\text{Type3}}(r;\epsilon,\sigma) & \text{if } \epsilon > 0 \\
     LJ_{\text{Type3}}(r;|\epsilon|,\sigma) & \text{if } \epsilon < 0
    \end{cases}

----

For the definition of the Lennard-Jones and Weeks-Chandler-Andersen potentials see the respective sections ( `LennardJones <LennardJones.html>`_ and `WCA <WCA.html>`_).

For all types:

* :math:`\epsilon` is the depth of the potential well (can be positive or negative)
* :math:`\sigma` is the distance at which the potential is zero
* :math:`r` is the distance between the particles

----

* **type**: ``NonBonded``, ``GeneralLennardJonesType1`` (or ``GeneralLennardJonesType2`` or ``GeneralLennardJonesType3``)
* **parameters**:
  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma

* **data**:

  * ``name_i``: ``string``: Type of the particle i
  * ``name_j``: ``string``: Type of the particle j
  * ``epsilon``: ``real``: Potential well depth :math:`[energy]`
  * ``sigma`` : ``real``: Zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "generalLJ":{
     "type":["NonBonded","GeneralLennardJonesType1"],
     "parameters":{
       "cutOffFactor":2.5,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", -1.2, 0.9],
       ["B", "B", 0.8, 1.1]
     ]
   }

.. tip::
    The GeneralLennardJones potential behaves like a standard Lennard-Jones potential for negative epsilon values and like a WCA potential for positive epsilon values. This allows for modeling both attractive and purely repulsive interactions within the same potential form.

.. note::
   For an explanation of how select the neighbors list (and the use of the `condition` parameter), see the common documentation for the non-bonded potentials.
