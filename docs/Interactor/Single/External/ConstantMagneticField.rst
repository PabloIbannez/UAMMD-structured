ConstantMagneticField
---------------------

The ConstantMagneticField potential models a uniform, constant magnetic field applied to magnetic particles in the system.

The magnetic field is given by:

.. math::

   \mathbf{B} = B_0 \vec{n}

where :math:`B_0` is the magnitude and :math:`\vec{n}` is the direction vector.

The energy of a magnetic dipole :math:`\mathbf{m}` in this field is:

.. math::

   U = -\mathbf{m} \cdot \mathbf{B}

And the torque on the dipole is:

.. math::

   \mathbf{\tau} = \mathbf{m} \times \mathbf{B}


----

* **type**: ``External``, ``ConstantMagneticField``
* **parameters**:

  * ``b0``: ``real``: Magnitude of the magnetic field :math:`[magnetic field]`, B_0
  * ``direction``: ``real3``: Direction vector of the magnetic field

Example:

.. code-block::

   "constantMagneticField":{
     "type":["External","ConstantMagneticField"],
     "parameters":{
       "b0": 1.0,
       "direction": [0.0, 0.0, 1.0]
     }
   }

.. note::
   This potential requires magnetic moment information for each particle. Ensure that the magnetic moments are properly set for each particle before using this potential.

.. tip::
   The ConstantMagneticField potential can be used to study the alignment and dynamics of magnetic particles or to model systems in strong external magnetic fields.
