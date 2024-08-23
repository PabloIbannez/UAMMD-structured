ExternalTabulated
-----------------

The ExternalTabulated potential allows for the application of an arbitrary external potential defined on a 3D grid. 
This is particularly useful for complex external fields that cannot be easily described analytically.

The potential is defined on a 3D grid with dimensions `nx`, `ny`, and `nz`. 
The simulation box (provided by the ensemble) is divided into `nx` x `ny` x `nz` grid points, and the potential is defined at each grid point.
The size of each grid cell the simulation box is divided into is given by:

.. math::
   \Delta x = \frac{L_x}{nx}, \quad \Delta y = \frac{L_y}{ny}, \quad \Delta z = \frac{L_z}{nz}

where `L_x`, `L_y`, and `L_z` are the lengths of the simulation box in the x, y, and z directions.

.. warning::
    The current version of the ExternalTabulated potential requires all the :math:`\Delta x`, :math:`\Delta y`, and :math:`\Delta z` to be **equal**.
    Otherwise, an error will be thrown.

An scaling factor `scale` can be provided. If so, the energy and force values are multiplied by this factor.

.. warning::
   The force values are interpolated from the given values for the force, not calculated from the energy values. (The derivative of the energy values is not used to calculate the force values.)
   Ensure that the force values are consistent with the energy values.

----

* **type**: ``External``, ``ExternalTabulated``
* **parameters**:

  * ``nx``, ``ny``, ``nz``: ``int``: Number of grid points in x, y, and z directions
  * ``scale``: ``real``: Scaling factor for energy and force (optional, default: 1.0)

* **data**:

  * ``i``, ``j``, ``k``: ``int``: Grid point indices
  * ``energy``: ``real``: Energy value at the grid point
  * ``force``: ``real3``: Force vector at the grid point

Example:

.. code-block::

   "externalTabulated":{
     "type":["External","ExternalTabulated"],
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

The potential and forces are interpolated from the provided grid values to particle positions using a 3-point interpolation scheme.

.. note::
   This potential is highly flexible and can represent complex external fields, but requires pre-computation of the energy and force values on a grid.

.. note::
   ExternalTabulated uses the interpolation provided by the Immerse Boundary Method (IBM) of UAMMD (Developed by Raul. P. Pelaez).
   Documentation for the IBM can be found at https://uammd.readthedocs.io/en/latest/ImmersedBoundary.html


.. warning::
   The tabulated data should cover the entire simulation box with sufficient resolution for accurate interpolation. Ensure that the grid spacing is fine enough to capture the relevant features of your external potential.

.. tip::
   The ``scale`` parameter can be used to easily adjust the strength of the external potential without recalculating the entire grid.
