SortStep
--------

The SortStep is a utility simulation step that sorts the particles in the system. This can be useful for optimizing memory access patterns and improving simulation performance.

----

* **type**: ``UtilsStep``, ``SortStep``
* **parameters**: None

Example:

.. code-block::

   "sortParticles":{
     "type":["UtilsStep","SortStep"],
     "parameters":{}
   }

.. note::
   Sorting particles can improve cache coherence and potentially speed up other operations in the simulation.

.. warning::
   Frequent sorting may have a performance cost, so use this step judiciously.
