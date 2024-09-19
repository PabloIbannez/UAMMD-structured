CompressibleInertialCoupling
----------------------------

The CompressibleInertialCoupling integrator implements a fluctuating hydrodynamics method for simulating particles coupled to a compressible fluid. This integrator is a wrapper around the UAMMD ICM_Compressible integrator.

For more details on the underlying method, please refer to the `UAMMD ICM_Compressible documentation <https://uammd.readthedocs.io/en/latest/Integrator/FluctuatingHydrodynamics.html#icm-compressible>`_.

----

* **type**: ``FluctuatingHydrodynamics``, ``CompressibleInertialCoupling``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``shearViscosity``: ``real``: Shear viscosity of the fluid :math:`[pressure \cdot time]`
  * ``bulkViscosity``: ``real``: Bulk viscosity of the fluid :math:`[pressure \cdot time]`
  * ``speedOfSound``: ``real``: Speed of sound in the fluid :math:`[distance/time]`
  * ``density``: ``real``: Density of the fluid :math:`[mass/distance^3]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]`
  * ``cellDimX``, ``cellDimY``, ``cellDimZ``: ``int``: Number of cells in each direction
  * ``initialCondition``: ``string``: Initial condition for the fluid (default: "none")

Example:

.. code-block::

   "compressibleICM":{
     "type":["FluctuatingHydrodynamics","CompressibleInertialCoupling"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "shearViscosity": 1.0,
       "bulkViscosity": 1.0,
       "speedOfSound": 10.0,
       "density": 1.0,
       "hydrodynamicRadius": 0.5,
       "cellDimX": 32,
       "cellDimY": 32,
       "cellDimZ": 32,
       "initialCondition": "none"
     }
   }

.. note::
   This integrator requires that the particle group contains all particles in the system.

.. warning::
   The ``initialCondition`` parameter currently only supports the "none" option, which initializes the fluid with uniform density and zero velocity.
