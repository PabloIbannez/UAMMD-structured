InfoStep
--------

The InfoStep is a utility simulation step that prints information about the simulation progress, including estimated time remaining and performance metrics.

The output is displayed through the standard output stream, and includes the following information:

- The current step
- ETA (estimated time of arrival) for the simulation to complete
- Mean FPS, mean number of steps completed per second

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
