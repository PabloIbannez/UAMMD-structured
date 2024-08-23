PolarizationTabulated
---------------------

The PolarizationTabulated potential is an extension of the `ExternalTabulated potential <ExternalTabulated.html>`_
that takes into account the polarizability of particles. 
It modifies the external field effect based on each particle's polarizability. In particular, the energy and force 
values are scaled by the polarizability of the particle (which is given by the type of the particle or defined in *State*).

----

* **type**: ``External``, ``PolarizationTabulated``
* **parameters**:

  * ``nx``, ``ny``, ``nz``: ``int``: Number of grid points in x, y, and z directions
  * ``scale``: ``real``: Scaling factor for energy and force (optional, default: 1.0)

* **data**:

  * ``i``, ``j``, ``k``: ``int``: Grid point indices
  * ``energy``: ``real``: Energy value at the grid point
  * ``force``: ``real3``: Force vector at the grid point

Example:

.. code-block::

   "polarizationTabulated":{
     "type":["External","PolarizationTabulated"],
     "parameters":{
       "nx":64, "ny":64, "nz":64,
       "scale":1.0
     },
     "labels":["i", "j", "k", "energy", "force"],
     "data":[
       [0, 0, 0, 1.0, [0.1, 0.0, 0.0]],
       [0, 0, 1, 0.9, [0.09, 0.0, 0.01]],
       ...
     ]
   }

.. note::
   This potential requires polarizability information for each particle. Ensure that the polarizability is properly set for each particle before using this potential.

.. warning::
   The tabulated data should cover the entire simulation box with sufficient resolution for accurate interpolation.
