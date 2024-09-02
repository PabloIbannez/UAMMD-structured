KimHummer
---------

The KimHummer potential is a coarse-grained model for protein-protein interactions developed by Kim and Hummer [Hummer2008]_.
It combines electrostatic interactions with a modified Lennard-Jones potential and includes solvent-accessible surface area (SASA) weighting.

The potential consists of two main terms:

.. math::

    U_{KH} = U_{el} + U_{np}

where:

* :math:`U_{el}` is the electrostatic term (Debye-Hückel potential)
* :math:`U_{np}` is the non-polar term (modified Lennard-Jones potential)

The electrostatic term is given by the Debye-Hückel potential:

.. math::

    U_{el} = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r} \frac{e^{-r/\lambda_D}}{r}

The non-polar term is a modified Lennard-Jones potential:

.. math::

    U_{np} = \epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} - 2\left(\frac{\sigma_{ij}}{r}\right)^6 \right]

Both terms are weighted by SASA-dependent factors for each residue.

----

* **type**: ``NonBonded``, ``KimHummer``
* **parameters**:

  * ``cutOffNPFactor``: ``real``: Interaction range for non-polar term as a multiple of sigma
  * ``cutOffDHFactor``: ``real``: Interaction range for electrostatics as a multiple of the Debye length
  * ``dielectricConstant``: ``real``: Relative permittivity of the medium :math:`\epsilon_r`
  * ``debyeLength``: ``real``: Debye length :math:`\lambda_D` :math:`[distance]`
  * ``sasaModel``: ``string``: SASA weighting model (options: "A", "B", "C", "D", "E", "F")
  * ``zeroEnergy``: ``real``: Energy shift parameter for non-polar term (default: 0.01) :math:`[energy]`

* **data**:

  * ``name_i``: ``string``: Type of residue i
  * ``name_j``: ``string``: Type of residue j
  * ``epsilon``: ``real``: LJ potential well depth :math:`[energy]`
  * ``sigma``: ``real``: LJ zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "kimHummer":{
     "type":["NonBonded","KimHummer"],
     "parameters":{
       "cutOffNPFactor":2.5,
       "cutOffDHFactor":3.0,
       "dielectricConstant":78.5,
       "debyeLength":1.0,
       "sasaModel":"B",
       "zeroEnergy":0.01
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["ALA", "ALA", -0.137, 5.0],
       ["ALA", "GLY", -0.068, 4.9],
       ["GLY", "GLY", 0.0, 4.8]
     ]
   }

.. note::
   The KimHummer potential requires charge, radius, and SASA information for each particle. Ensure these properties are properly set before using this potential.

.. [Hummer2008] Kim, Y. C., & Hummer, G. (2008). Coarse-grained models for simulations of multiprotein complexes: application to ubiquitin binding. Journal of molecular biology, 375(5), 1416-1433.
