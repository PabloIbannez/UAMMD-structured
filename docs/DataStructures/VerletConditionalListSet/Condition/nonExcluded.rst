non excluded
------------

.. note::

   For a description of the common aspects of conditions, including their usage and shared parameters (``cutOff``, ``cutOffVerletFactor``), please refer to the `general conditions documentation <index.html>`_.

This condition generates one Verlet list only.

It is ensured that if a pairs is given in the exclusions list, it will not be in the Verlet list.
To access the Verlet list, use the label ``"nonExcluded"`` in the ``"condition"`` option of the potential.

To introduce the exclusion list, we use the ``"labels"`` and ``"data"`` fields of the **DataEntry**.
In the labels field, we have to specify two columns, named ``"id"`` and ``"id_list"``.
The ``"id"`` column contains the id of the particle, and the ``"id_list"`` 
column contains the list of ids of the particles that are excluded from the Verlet list.

For example if we want the particle with ``id 0`` do not interact with the particles with ``id 1,2,3`` 
we have to specify the following data:

.. code-block::

   "data":[
      [0,[1,2,3]],
      [1,[0]],
      [2,[0]],
      [3,[0]],
      ["...","..."]
   ]

.. warning::

   If the particle with ``id 0`` is excluded from the Verlet list of the particle with ``id 1``, 
   it does not mean that the particle with ``id 1`` is excluded from the Verlet list of the particle with ``id 0``.
   This has to be written explicitly in the data. If the given exclusion list is not symmetric, the code will raise an error.


----

* **type**: ``VerletConditionalListSet``, ``nonExcluded``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**:

   .. list-table::
      :widths: 25 25
      :header-rows: 1
      :align: center

      * - id
        - id_list
      * - ``int``
        - [``int``, ``int`` , ...]

----

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","nonExcluded"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
     "labels":["id","id_list"],
     "data":[
        [0,[1,2,3]],
        [1,[0,2,3]],
        [2,[0,1,3]],
        [3,[0,1,2]],
        ["...","..."]
     ]
   }
