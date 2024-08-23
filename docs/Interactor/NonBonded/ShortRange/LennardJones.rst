LennardJones
------------

The Lennard-Jones potential models both attractive and repulsive forces between particles. 
It's commonly used for modeling non specific interactions between particles.

There are three types of Lennard-Jones potentials implemented:

LennardJonesType1
~~~~~~~~~~~~~~~~~

.. math::

    U = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]

LennardJonesType2
~~~~~~~~~~~~~~~~~

.. math::

    U = \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 2\left(\frac{\sigma}{r}\right)^6 \right]

LennardJonesType3
~~~~~~~~~~~~~~~~~

.. math::

    U = \epsilon \left[ 5\left(\frac{\sigma}{r}\right)^{12} - 6\left(\frac{\sigma}{r}\right)^{10} \right]

For all types:

* :math:`\epsilon` is the depth of the potential well
* :math:`\sigma` is the distance at which the potential is zero
* :math:`r` is the distance between the particles

----

* **type**: ``NonBonded``, ``LennardJonesType1`` (or ``LennardJonesType2`` or ``LennardJonesType3``)

* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma

* **data**:

  * ``name_i``: ``string``: Type of the particle i
  * ``name_j``: ``string``: Type of the particle j
  * ``epsilon``: ``real``: Potential well depth :math:`[energy]`
  * ``sigma``  : ``real``: Zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "lj":{
     "type":["NonBonded","LennardJonesType1"],
     "parameters":{
       "cutOffFactor":2.5,
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
   For an explanation of how select the neighbors list (and the use of the `condition` parameter), see the common documentation for the non-bonded potentials.
