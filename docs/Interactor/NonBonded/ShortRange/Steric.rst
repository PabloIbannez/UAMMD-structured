Steric
------

The Steric potential is used to model repulsive interactions between particles, typically representing excluded volume effects. It is implemented in two forms: Steric6 and Steric12.

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
* :math:`r` is the distance between particles

----

* **type**: ``NonBonded``, ``Steric6`` or ``Steric12``
* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Interaction strength :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`

Example:

.. code-block::

   "steric":{
     "type":["NonBonded","Steric6"],
     "parameters":{
       "cutOffFactor":2.5,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 0.8, 0.9],
       ["B", "B", 1.2, 1.1]
     ]
   }

.. tip::
   The Steric potential is purely repulsive and is often used to model the core repulsion between particles or to prevent overlap in coarse-grained simulations.
