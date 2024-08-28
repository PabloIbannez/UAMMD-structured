AFMMeasure
----------

The AFMMeasure step simulates an Atomic Force Microscopy (AFM) experiment, recording force-distance curves and other relevant data.

The output file format depends on the ``mode`` parameter:

**Standard mode**:

.. code-block::

   #Step                    X                    F
   [step number]     [tip height]     [tip deflection force]

Where:

- ``Step`` is the simulation step number
- ``X`` is the tip height (tip position minus tip radius)
- ``F`` is the tip deflection force

**Verbose mode**:

.. code-block::

   #Step        ChipPos        TipPos    tipDeflection    tipForce    sampleForce    tipDeflectionForce
   [step]       [chip z]       [tip z]    [deflection]    [tip Fz]    [sample Fz]    [deflection F]

Where:

- ``Step`` is the simulation step number
- ``ChipPos`` is the z-position of the AFM chip
- ``TipPos`` is the z-position of the AFM tip
- ``tipDeflection`` is the deflection of the tip from the chip position
- ``tipForce`` is the z-component of the force on the tip
- ``sampleForce`` is the z-component of the force on the sample
- ``tipDeflectionForce`` is the force due to tip deflection

All values are written in scientific notation with 24 character width for alignment.

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
       "intervalStep": 10000,
       "outputFilePath": "afm_data.dat",
       "mode": "Verbose",
       "tipType": "TIP"
     }
   }

.. note::
   This step requires an AFM interactor to be present in the simulation.

.. tip::
   The "Verbose" mode provides detailed output including chip position, tip position, and forces, while "Standard" mode gives a more concise force-distance output.
