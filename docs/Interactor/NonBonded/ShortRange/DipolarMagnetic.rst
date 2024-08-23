DipolarMagnetic
---------------

The DipolarMagnetic potential models the interaction between magnetic dipoles. It is useful for simulating systems of magnetic nanoparticles or other materials with significant magnetic moments.

The potential energy between two magnetic dipoles is given by:

.. math::

    U = -\frac{\mu_0}{4\pi r^3} [3(\mathbf{m}_i \cdot \hat{\mathbf{r}})(\mathbf{m}_j \cdot \hat{\mathbf{r}}) - \mathbf{m}_i \cdot \mathbf{m}_j]

where:

* :math:`\mu_0` is the vacuum permeability
* :math:`r` is the distance between dipoles
* :math:`\mathbf{m}_i, \mathbf{m}_j` are the magnetic moments of particles i and j
* :math:`\hat{\mathbf{r}}` is the unit vector pointing from particle i to j

----

* **type**: ``NonBonded``, ``DipolarMagnetic``
* **parameters**:

  * ``cutOff``: ``real``: Interaction cutoff distance :math:`[distance]`
  * ``permeability``: ``real``: Magnetic permeability of the medium

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j

Example:

.. code-block::

   "dipolarMagnetic":{
     "type":["NonBonded","DipolarMagnetic"],
     "parameters":{
       "cutOff":10.0,
       "permeability":1.25663706e-6,
       "condition":"all"
     },
     "labels":["name_i", "name_j"],
     "data":[
       ["A", "A"],
       ["A", "B"],
       ["B", "B"]
     ]
   }

.. note::
   The DipolarMagnetic potential requires magnetic moment information for each particle. Ensure that the magnetic moments are properly set for each particle before using this potential.

.. warning::
   Long-range corrections may be necessary for accurate results in systems with significant long-range dipolar interactions.
