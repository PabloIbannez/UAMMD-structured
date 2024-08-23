all
----

.. note::

   For a description of the common aspects of conditions, including their usage and shared parameters (``cutOff``, ``cutOffVerletFactor``), please refer to the `general conditions documentation <index.html>`_.


Using this condition we create a Verlet list with defines one condition only, "all".
With the "all" condition all possible pairs of particles are considered for the Verlet list (using this condition is equivalent to create a regular Verlet list).

----

* **type**: ``VerletConditionalListSet``, ``all``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.
  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

----

Example:

.. code-block:: 

   "entryName":{
     "type":["VerletConditionalListSet","all"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
   }

