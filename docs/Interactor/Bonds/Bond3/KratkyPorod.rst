KratkyPorod
-----------

The Kratky-Porod potential, also known as the worm-like chain model, is used to describe the bending energy of semi-flexible polymers. It penalizes deviations from a straight chain configuration.

.. math::

    U = K(1 + \cos\theta)

where:

* :math:`K` is the bending stiffness
* :math:`\theta` is the angle between two consecutive segments

----

* **type**: ``Bond3``, ``KratkyPorod``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle
  * ``K``: ``real``: Bending stiffness :math:`[energy]`

Example:

.. code-block::

   "kratkyPorodBonds":{
     "type":["Bond3","KratkyPorod"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k", "K"],
     "data":[[0, 1, 2, 10.0],
             [1, 2, 3, 10.0],
             [2, 3, 4, 10.0]]
   }

KratkyPorodCommon_K
~~~~~~~~~~~~~~~~~~~

KratkyPorod bonds variant with a common bending stiffness (``K``) for all bonds.

----

* **type**: ``Bond3``, ``KratkyPorodCommon_K``
* **parameters**:

  * ``K``: ``real``: Common bending stiffness for all bonds :math:`[energy]`

* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle

Example:

.. code-block::

   "kratkyPorodBondsCommonK":{
     "type":["Bond3","KratkyPorodCommon_K"],
     "parameters":{"K":10.0},
     "labels":["id_i", "id_j", "id_k"],
     "data":[[0, 1, 2],
             [1, 2, 3],
             [2, 3, 4]]
   }
