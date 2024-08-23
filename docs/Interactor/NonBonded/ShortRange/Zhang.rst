Zhang
-----

The Zhang potential is a specialized interaction model developed for coarse-grained simulations of anisotropic particles, particularly useful for modeling protein-protein interactions with orientation dependence.

The potential combines a distance-dependent term with an orientation-dependent term:

.. math::

    U = U_{\text{dist}}(r) \cdot U_{\text{orient}}(\mathbf{r}, \mathbf{n}_i, \mathbf{n}_j)

where:

* :math:`r` is the distance between particles
* :math:`\mathbf{n}_i, \mathbf{n}_j` are the orientation vectors of particles i and j

The distance-dependent term :math:`U_{\text{dist}}(r)` is typically a combination of attractive and repulsive parts, while the orientation-dependent term :math:`U_{\text{orient}}(\mathbf{r}, \mathbf{n}_i, \mathbf{n}_j)` modulates the interaction based on the relative orientations of the particles.

----

* **type**: ``NonBonded``, ``Zhang``
* **parameters**:

  * ``cutOffNPFactor``: ``real``: Interaction range as a multiple of sigma
  * ``cutOffDHFactor``: ``real``: Interaction range for orientation-dependent term
  * ``axis``: ``real3``: Reference axis for orientation (default: {0, 0, 1})

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Strength of interaction :math:`[energy]`
  * ``sigma``: ``real``: Characteristic distance :math:`[distance]`
  * ``mu``: ``real``: Anisotropy parameter
  * ``theta``: ``real``: Preferred angle :math:`[angle]`

Example:

.. code-block::

   "zhang":{
     "type":["NonBonded","Zhang"],
     "parameters":{
       "cutOffNPFactor":2.5,
       "cutOffDHFactor":3.0,
       "axis":[0, 0, 1],
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma", "mu", "theta"],
     "data":[
       ["A", "A", 1.0, 1.0, 0.5, 1.57],
       ["A", "B", 1.2, 0.9, 0.6, 1.2],
       ["B", "B", 0.8, 1.1, 0.4, 1.8]
     ]
   }

.. note::
   The Zhang potential requires orientation information for each particle. Ensure that the orientation (director) is properly set for each particle before using this potential.

.. warning::
   The Zhang potential is specifically designed for anisotropic interactions and may not be suitable for isotropic systems without modification.
