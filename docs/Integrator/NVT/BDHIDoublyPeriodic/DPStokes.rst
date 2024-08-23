DPStokes
--------

The DPStokes integrator implements Brownian Dynamics with Hydrodynamic Interactions for doubly periodic systems using the Doubly Periodic Stokes (DPStokes) method. This integrator is a wrapper around the UAMMD DPStokesSlab_ns::DPStokesIntegrator.

For more details on the underlying method, please refer to the `UAMMD DPStokes documentation <https://uammd.readthedocs.io/en/latest/Integrators.html#dpstokes>`_.

----

* **type**: ``BDHIDoublyPeriodic``, ``DPStokes``
* **parameters**:

  * ``timeStep``: ``real``: Time step :math:`[time]`
  * ``temperature``: ``real``: Temperature of the system :math:`[energy]`
  * ``viscosity``: ``real``: Viscosity of the fluid :math:`[mass/(distance \cdot time)]`
  * ``hydrodynamicRadius``: ``real``: Hydrodynamic radius of the particles :math:`[distance]`
  * ``nx``, ``ny``: ``int``: Number of grid points in x and y directions
  * ``nz``: ``int``: Number of grid points in z direction (optional)
  * ``Lx``, ``Ly``: ``real``: Box size in x and y directions :math:`[distance]`
  * ``H``: ``real``: Height of the simulation box :math:`[distance]`
  * ``w``: ``real``: Width of the Gaussian kernel :math:`[distance]`
  * ``beta``: ``real``: Regularization parameter for the Gaussian kernel (optional)
  * ``alpha``: ``real``: Regularization parameter for the wall (optional)
  * ``w_d``: ``real``: Width of the Gaussian kernel for the derivative (optional)
  * ``beta_d``: ``real``: Regularization parameter for the derivative kernel (optional)
  * ``alpha_d``: ``real``: Regularization parameter for the wall derivative (optional)
  * ``mode``: ``string``: Wall mode ("none", "bottom", or "slit")
  * ``tolerance``: ``real``: Tolerance for the iterative solver (default: 1e-7)

Example:

.. code-block::

   "dpstokes":{
     "type":["BDHIDoublyPeriodic","DPStokes"],
     "parameters":{
       "timeStep": 0.01,
       "temperature": 1.0,
       "viscosity": 1.0,
       "hydrodynamicRadius": 0.5,
       "nx": 64,
       "ny": 64,
       "Lx": 32.0,
       "Ly": 32.0,
       "H": 10.0,
       "w": 0.5,
       "mode": "slit",
       "tolerance": 1e-6
     }
   }

.. note::
   This integrator is suitable for systems with periodic boundary conditions in two dimensions and various boundary conditions in the third dimension.

.. warning::
   This integrator requires that the particle group contains all particles in the system.
