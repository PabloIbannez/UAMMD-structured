Clashed
-------

The Clashed potential is a soft repulsive interaction designed to prevent overlap between particles while allowing for some degree of interpenetration. It is particularly useful in coarse-grained simulations where hard-core repulsions may be too restrictive.

The potential energy is given by:

.. math::

    U = \begin{cases}
    \lambda (d^2 - r^2)^2 & \text{if } r < d \\
    0 & \text{if } r \geq d
    \end{cases}

where:

* :math:`\lambda` is the strength of the interaction
* :math:`d = \gamma(r_i + r_j)` is the interaction distance
* :math:`r_i, r_j` are the radii of particles i and j
* :math:`\gamma` is a scaling factor for the interaction distance
* :math:`r` is the distance between particle centers

----

* **type**: ``NonBonded``, ``Clashed``
* **parameters**:

  * ``lambda``: ``real``: Strength of the interaction :math:`[energy]/[distance]^4`
  * ``gamma``: ``real``: Scaling factor for the interaction distance

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j

Example:

.. code-block::

   "clashed":{
     "type":["NonBonded","Clashed"],
     "parameters":{
       "lambda":1.0,
       "gamma":1.1,
       "condition":"all"
     },
     "labels":["name_i", "name_j"],
     "data":[
       ["A", "A"],
       ["A", "B"],
       ["B", "B"]
     ]
   }

.. note::
   The Clashed potential uses the radius information stored in the particle data. Ensure that the radii are properly set for each particle before using this potential.

.. warning::
   While the Clashed potential allows for some overlap between particles, it may not be suitable for systems where maintaining strict excluded volume is critical.
