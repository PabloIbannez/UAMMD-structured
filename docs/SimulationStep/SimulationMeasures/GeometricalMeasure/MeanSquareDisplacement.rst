MeanSquareDisplacement
----------------------

The MeanSquareDisplacement step calculates and records the mean square displacement (MSD) of particles relative to their initial positions.

----

* **type**: ``GeometricalMeasure``, ``MeanSquareDisplacement``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

Example:

.. code-block::

   "msd":{
     "type":["GeometricalMeasure","MeanSquareDisplacement"],
     "parameters":{
       "outputFilePath": "msd.dat"
     }
   }

.. note::
   The output file will contain the step number and the mean square displacement value.

.. tip::
   The MSD is a key measure for studying diffusion and particle mobility in the system.
