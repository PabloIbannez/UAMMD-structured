DynamicallyBondedPatchyParticles
---------------------------------

This **Fundamental** is used by pathcy particles when the **Single Bond Per Patch** (SBPP)  is used.
This component defines the energy threshold that determines if a bond (between two patches) is formed or not.


* **type**: ``Fundamental``, ``DynamicallyBondedPatchyParticles``.
* **parameters**:

  * ``energyThreshold`` : ``float``, *optional*, default: 0.0.

* **data**: ``None``.

----

Example:

.. code-block:: json

   "entryName": {
     "type": ["Fundamental", "DynamicallyBondedPatchyParticles"],
     "parameters": {
       "energyThreshold": -1.0
     }
   }

.. note::
   Using this fundamental makes only sense whem we are using a patchy particle model with the SBPP model.
