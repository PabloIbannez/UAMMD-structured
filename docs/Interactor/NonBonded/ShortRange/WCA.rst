WCA (Weeks-Chandler-Andersen)
-----------------------------

The WCA potential is a shifted and truncated version of the Lennard-Jones potential, providing a purely repulsive interaction. It is often used to model excluded volume effects in soft matter systems.

There are three types of WCA potentials implemented:

WCAType1
~~~~~~~~

.. math::

    U = \begin{cases}
    4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] + \epsilon & \text{if } r < 2^{1/6}\sigma \\
    0 & \text{if } r \geq 2^{1/6}\sigma
    \end{cases}

WCAType2
~~~~~~~~

.. math::

    U = \begin{cases}
    \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 2\left(\frac{\sigma}{r}\right)^6 \right] + \epsilon & \text{if } r < \sigma \\
    0 & \text{if } r \geq \sigma
    \end{cases}

WCAType3
~~~~~~~~

.. math::

    U = \begin{cases}
    \epsilon \left[ 5\left(\frac{\sigma}{r}\right)^{12} - 6\left(\frac{\sigma}{r}\right)^{10} \right] + \epsilon & \text{if } r < \sigma \\
    0 & \text{if } r \geq \sigma
    \end{cases}

For all types:

* :math:`\epsilon` is the strength of the interaction
* :math:`\sigma` is the characteristic distance
* :math:`r` is the distance between particles

----

* **type**: ``NonBonded``, ``WCAType1`` (or ``WCAType2`` or ``WCAType3``)
* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma (typically 1.122462048309373 for Type1 and Type2, 1.0747892746492746 for Type3)

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Interaction strength :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`

Example:

.. code-block::

   "wca":{
     "type":["NonBonded","WCAType1"],
     "parameters":{
       "cutOffFactor":1.5,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 0.8, 0.9],
       ["B", "B", 1.2, 1.1]
     ]
   }

.. note::
   The WCA potential provides a continuous, purely repulsive interaction that is computationally efficient and widely used in molecular dynamics simulations of soft matter systems.
