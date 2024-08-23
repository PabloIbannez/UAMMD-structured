LennardJones
------------

The Lennard-Jones potential models both attractive and repulsive forces between particles. It's commonly used for van der Waals interactions but can also be applied as a bond potential for specific applications.

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

   U(r) = \epsilon \left[ 5\left(\frac{\sigma}{r}\right)^{12} - 6\left(\frac{\sigma}{r}\right)^{10} \right]

For all types:

* :math:`\epsilon` is the depth of the potential well
* :math:`\sigma` is the distance at which the potential is zero
* :math:`r` is the distance between the particles

----

* **type**: ``Bond2``, ``LennardJonesType1`` (or ``LennardJonesType2`` or ``LennardJonesType3``)
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``epsilon``: ``real``: Potential well depth :math:`[energy]`
  * ``sigma``  : ``real``: Zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "lennardJonesBonds":{
     "type":["Bond2","LennardJonesType1"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.0, 1.0],
             [1, 2, 1.2, 0.9],
             [2, 3, 0.8, 1.1]]
   }

LennardJonesType1Common_epsilon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LennardJones bonds variant with a common epsilon for all bonds. Similar variants exist for Type2 and Type3.

----

* **type**: ``Bond2``, ``LennardJonesType1Common_epsilon``
* **parameters**:

  * ``epsilon``: ``real``: Common potential well depth for all bonds :math:`[energy]`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``sigma``: ``real``: Zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "lennardJonesBondsCommonEpsilon":{
     "type":["Bond2","LennardJonesType1Common_epsilon"],
     "parameters":{"epsilon":1.0},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 1.0],
             [1, 2, 0.9],
             [2, 3, 1.1]]
   }
