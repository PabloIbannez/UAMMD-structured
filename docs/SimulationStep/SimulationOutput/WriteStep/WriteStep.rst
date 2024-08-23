WriteStep
---------

The WriteStep is a simulation step that writes the system state to a file at specified intervals during the simulation. It supports various output formats for particle positions, velocities, and other properties.

----

* **type**: ``WriteStep``, ``WriteStep``
* **parameters**:

  * ``outputFilePath``: ``string``: Path and base name for the output file
  * ``outputFormat``: ``string``: Format of the output file (e.g., "xyz", "pdb", "lammpstrj")
  * ``pbc``: ``bool``: Whether to apply periodic boundary conditions (optional, default: true)
  * ``append``: ``bool``: Whether to append to an existing file (optional, default: false)

Example:

.. code-block::

   "writeTrajectory":{
     "type":["WriteStep","WriteStep"],
     "parameters":{
       "outputFilePath": "trajectory",
       "outputFormat": "xyz",
       "pbc": true,
       "append": false
     }
   }

.. note::
   The actual file name will be {outputFilePath}.{outputFormat}

.. tip::
   Available output formats include "coord", "sp", "spo", "spf", "xyz", "pdb", "itpv", "itpd", "dcd", "lammpstrj", "vel", "magnet", "xyzm", "spm", "svv", "svvm", "svvma".
