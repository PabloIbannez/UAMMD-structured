WritePatchyParticlesStep
------------------------

The WritePatchyParticlesStep is an extension of the WriteStep that specifically handles systems with patchy particles. It writes both the core particles and their associated patches to the output file.

----

* **type**: ``WriteStep``, ``WritePatchyParticlesStep``
* **parameters**:

  * ``outputFilePath``: ``string``: Path and base name for the output file
  * ``outputFormat``: ``string``: Format of the output file (e.g., "xyz", "pdb", "lammpstrj")
  * ``pbc``: ``bool``: Whether to apply periodic boundary conditions (optional, default: true)
  * ``append``: ``bool``: Whether to append to an existing file (optional, default: false)

Example:

.. code-block::

   "writePatchyTrajectory":{
     "type":["WriteStep","WritePatchyParticlesStep"],
     "parameters":{
       "outputFilePath": "patchy_trajectory",
       "outputFormat": "xyz",
       "pbc": true,
       "append": false
     }
   }

.. note::
   This step automatically detects and processes PatchyParticles interactors in the system.

.. warning::
   Ensure that the PatchyParticles interactors are properly set up in your simulation before using this step.

.. tip::
   The WritePatchyParticlesStep is particularly useful for visualizing and analyzing systems with complex, patchy particle structures.
