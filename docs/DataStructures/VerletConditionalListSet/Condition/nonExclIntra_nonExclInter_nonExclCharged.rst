non excl intra - non excl inter - non excl charged
--------------------------------------------------

.. note::

   For a description of the common aspects of conditions, including their usage and shared parameters (``cutOff``, ``cutOffVerletFactor``), please refer to the `general conditions documentation <index.html>`_.

This condition creates three diffirent lists. All of them are non-excluded.
This means that is ensured if a pairs of particles is given in the exclusion list,
this pairs will not be in any of the three lists.

The three lists are:

1. A list with pairs of particles with the same ``modelId`` (and not excluded). This list is labeled as ``intra``.
2. A list with pairs of particles with different ``modelId`` (and not excluded). This list is labeled as ``inter``.
3. A list with pairs of particles where both particles have a charge (and not excluded). This list is labeled as ``charged``.

.. note::
   A detailed description of how to set the exclusion list can be found in the `documentation of nonExcluded condition <nonExcluded.html>`_.

----

* **type**: ``VerletConditionalListSet``, ``nonExclIntra_nonExclInter_nonExclCharged``.
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
     "type":["VerletConditionalListSet","nonExclIntra_nonExclInter_nonExclCharged"]
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
