ACMagneticField
---------------

The ACMagneticField potential models an alternating (AC) magnetic field applied to magnetic particles. It extends the ConstantMagneticField potential by adding time-dependent oscillation.

The magnetic field at time :math:`t` is given by:

.. math::

   \mathbf{B}(t) = B_0 \sin(2\pi f t + \phi) \hat{n}

where :math:`B_0` is the amplitude, :math:`f` is the frequency, :math:`\phi` is the phase, and :math:`\hat{n}` is the normalized direction vector.


----

* **type**: ``External``, ``ACMagneticField``
* **parameters**:

  * ``b0``: ``real``: Amplitude of the magnetic field :math:`[magnetic field]`
  * ``direction``: ``real3``: Direction vector of the magnetic field (will be normalized)
  * ``frequency``: ``real``: Frequency of the AC field :math:`[1/time]`
  * ``phase``: ``real``: Initial phase of the AC field :math:`[angle]` (optional, default: 0)

Example:

.. code-block::

   "acMagneticField":{
     "type":["External","ACMagneticField"],
     "parameters":{
       "b0": 1.0,
       "direction": [0.0, 0.0, 1.0],
       "frequency": 1.0,
       "phase": 0.0
     }
   }

.. note::
   This potential requires magnetic moment information for each particle. Ensure that the magnetic moments are properly set for each particle before using this potential.

.. warning::
   The ACMagneticField potential may require small time steps to accurately resolve the oscillations, especially at high frequencies.
