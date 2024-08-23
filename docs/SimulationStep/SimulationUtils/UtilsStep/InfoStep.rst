InfoStep
--------

The InfoStep is a utility simulation step that prints information about the simulation progress, including estimated time remaining and performance metrics.

----

* **type**: ``UtilsStep``, ``InfoStep``
* **parameters**: None

Example:

.. code-block::

   "printInfo":{
     "type":["UtilsStep","InfoStep"],
     "parameters":{}
   }

.. note::
   This step is useful for monitoring long-running simulations and estimating completion times.

.. tip::
   The InfoStep can help identify performance bottlenecks or unexpected slowdowns during the simulation.
