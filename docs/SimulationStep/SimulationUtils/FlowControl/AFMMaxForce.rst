AFMMaxForce
-----------

The AFMMaxForce is a flow control simulation step designed for Atomic Force Microscopy (AFM) simulations. It stops the simulation when all AFM tips have reached a specified maximum force.

----

* **type**: ``FlowControl``, ``AFMMaxForce``
* **parameters**:

  * ``maxForce``: ``real``: Maximum force threshold for AFM tips :math:`[force]`
  * ``tipType``: ``string``: Type name for AFM tip particles (default: "TIP")

Example:

.. code-block::

   "afmMaxForce":{
     "type":["FlowControl","AFMMaxForce"],
     "parameters":{
       "maxForce": 100.0,
       "tipType": "TIP"
     }
   }

.. note::
   This step automatically detects and processes AFM interactors in the system.

.. warning::
   Ensure that AFM interactors are properly set up in your simulation before using this step.
