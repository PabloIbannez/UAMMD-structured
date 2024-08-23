AFMMeasure
----------

The AFMMeasure step simulates an Atomic Force Microscopy (AFM) experiment, recording force-distance curves and other relevant data.

----

* **type**: ``ExperimentMeasures``, ``AFMMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file
  * ``mode``: ``string``: Output mode, either "Verbose" or "Standard" (default: "Standard")
  * ``tipType``: ``string``: Type name for AFM tip particles (default: "TIP")

Example:

.. code-block::

   "afmMeasure":{
     "type":["ExperimentMeasures","AFMMeasure"],
     "parameters":{
       "outputFilePath": "afm_data.dat",
       "mode": "Verbose",
       "tipType": "TIP"
     }
   }

.. note::
   This step requires an AFM interactor to be present in the simulation.

.. tip::
   The "Verbose" mode provides detailed output including chip position, tip position, and forces, while "Standard" mode gives a more concise force-distance output.
