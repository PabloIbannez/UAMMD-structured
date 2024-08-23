intramolecular-intermolecular
-----------------------------

.. note::

   For a description of the common aspects of conditions, including their usage and shared parameters (``cutOff``, ``cutOffVerletFactor``), please refer to the `general conditions documentation <index.html>`_.

Using this condition two Verlet list are generated.

1. The first list contains pairs with the same ``modelId``, this list is labeled as ``"intra"``.
2. The second list contains pairs with different ``modelId``, this list is labeled as ``"inter"``.

----

* **type**: ``VerletConditionalListSet``, ``intra_inter``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

----

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","intra_inter"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
   }
