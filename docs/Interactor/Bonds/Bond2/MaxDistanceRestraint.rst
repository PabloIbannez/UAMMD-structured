MaxDistanceRestraint
--------------------

The MaxDistanceRestraint is a bond potential that applies a harmonic restraint when the distance between two particles exceeds a specified maximum distance. This can be useful for maintaining the overall structure of a system while allowing some flexibility.

.. math::

    U = \begin{cases}
    0 & \text{if } r \leq r_{\text{max}} \\
    \frac{1}{2}K(r - r_{\text{max}})^2 & \text{if } r > r_{\text{max}}
    \end{cases}

where:

* :math:`K` is the spring constant
* :math:`r` is the current distance between the bonded particles
* :math:`r_{\text{max}}` is the maximum allowed distance

----

* **type**: ``Bond2``, ``MaxDistanceRestraint``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``K``: ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``maxDistance``: ``real``: Maximum allowed distance :math:`[distance]`

Example:

.. code-block::

   "maxDistanceRestraints":{
     "type":["Bond2","MaxDistanceRestraint"],
     "parameters":{},
     "labels":["id_i", "id_j", "K", "maxDistance"],
     "data":[[0, 1, 100.0, 2.0],
             [1, 2, 100.0, 2.5],
             [2, 3, 100.0, 2.2]]
   }
