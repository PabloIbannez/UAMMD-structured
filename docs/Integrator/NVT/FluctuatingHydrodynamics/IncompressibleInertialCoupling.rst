IncompressibleInertialCoupling
------------------------------

The IncompressibleInertialCoupling integrator implements a fluctuating hydrodynamics method for simulating particles coupled to an incompressible fluid. This integrator is a wrapper around the UAMMD ICM integrator.

For more details on the underlying method, please refer to the `UAMMD ICM documentation <https://uammd.readthedocs.io/en/latest/Integrators.html#icm>`_.

----

* **type**: ``FluctuatingHydrodynamics``, ``IncompressibleInertialCoupling``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[pressure \cdot time]`
  * ``density``: ``real``: Density of the fluid :math:`[mass/distance^3]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]`
  * ``cellDimX``, ``cellDimY``, ``cellDimZ``: ``int``: Number of cells in each direction (optional)
  * ``sumThermalDrift``: ``bool``: Whether to sum the thermal drift term (default: false)
  * ``removeTotalMomentum``: ``bool``: Whether to remove the total momentum of the system (default: true)

Example:

.. code-block::

   "incompressibleICM":{
     "type":["FluctuatingHydrodynamics","IncompressibleInertialCoupling"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "density": 1.0,
       "hydrodynamicRadius": 0.5,
       "cellDimX": 32,
       "cellDimY": 32,
       "cellDimZ": 32,
       "sumThermalDrift": false,
       "removeTotalMomentum": true
     }
   }

.. note::
   If ``cellDimX``, ``cellDimY``, and ``cellDimZ`` are not provided, the integrator will automatically determine the cell dimensions.

.. warning::
   This integrator requires that the particle group contains all particles in the system.
