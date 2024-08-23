Electrostatics
--------------

The Electrostatics interactor is a wrapper around the UAMMD SpectralEwaldPoisson module, implementing long-range electrostatic interactions using the Spectral Ewald Poisson method. This interactor is designed to handle charged particles in periodic systems efficiently.

For more details on the underlying method, please refer to the `UAMMD SpectralEwaldPoisson documentation <https://uammd.readthedocs.io/en/latest/Interactors.html#spectralewaldpoisson>`_.

----

* **type**: ``LongRange``, ``Electrostatics``
* **parameters**:

  * ``dielectricConstant``: ``real``: Dielectric constant of the medium
  * ``gaussianWidth``: ``real``: Width of the Gaussian used in the Spectral Ewald method :math:`[distance]`
  * ``tolerance``: ``real``: Tolerance for the Spectral Ewald method (default: 1e-4)
  * ``splitThreshold``: ``real``: Threshold for splitting the computation between real and Fourier space (default: 0.0, automatic)

Example:

.. code-block::

   "electrostatics":{
     "type":["LongRange","Electrostatics"],
     "parameters":{
       "dielectricConstant": 78.5,
       "gaussianWidth": 1.0,
       "tolerance": 1e-5,
       "splitThreshold": 0.5
     }
   }

.. note::
   This interactor uses the UAMMD SpectralEwaldPoisson module to compute long-range electrostatic interactions efficiently in periodic systems.

.. warning::
   The Electrostatics interactor currently only supports energy and force calculations. Other computables such as stress, hessian, and magnetic field are not implemented and will generate warnings if requested.

Functionality:
   - Computes electrostatic energies and forces for charged particles in periodic systems.
   - Uses the Spectral Ewald Poisson method for efficient long-range interactions.
   - Automatically adjusts to changes in the system's dielectric constant.

Limitations:
   - The simulation box cannot be changed after the interactor is initialized.
   - Stress, hessian, lambda derivative, magnetic field, and pairwise force calculations are not supported.

The Electrostatics interactor is particularly useful for simulations of charged particles in periodic systems, such as electrolytes, ionic solutions, or charged colloids. It provides an efficient method for handling long-range electrostatic interactions, which is crucial for accurate simulations of these systems.
