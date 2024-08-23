intermolecular different type
-----------------------------

.. note::

   For a description of the common aspects of conditions, including their usage and shared parameters (``cutOff``, ``cutOffVerletFactor``), please refer to the `general conditions documentation <index.html>`_.

This condition creates one Verlet list only. This Verlet list can be
requested using the identifier ``"inter"`` in the "condition" option of the 
potential.

All the pairs contains only particles with different ``modelId`` values (intermolecular) and
with different ``type`` (different type).

----

* **type**: ``VerletConditionalListSet``, ``interDifferentType``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.
  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

----

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","interDifferentType"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     }
   }
