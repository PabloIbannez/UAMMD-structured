KimHummer
---------

The KimHummer potential is a specialized interaction model developed by Kim and Hummer for coarse-grained simulations of protein-protein interactions. It combines electrostatic interactions with a modified Lennard-Jones potential and includes a term based on the solvent-accessible surface area (SASA).

The potential is composed of three terms:

.. math::

    U_{KH} = U_{el} + U_{LJ} + U_{SASA}

where:

* :math:`U_{el}` is the electrostatic term (Debye-Hückel potential)
* :math:`U_{LJ}` is a modified Lennard-Jones term
* :math:`U_{SASA}` is the SASA-dependent term

The electrostatic term is given by the Debye-Hückel potential:

.. math::

    U_{el} = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r} \frac{e^{-r/\lambda_D}}{r}

The modified Lennard-Jones term is:

.. math::

    U_{LJ} = \epsilon_{ij} \left[ 13\left(\frac{\sigma_{ij}}{r}\right)^{12} - 18\left(\frac{\sigma_{ij}}{r}\right)^{10} + 4\left(\frac{\sigma_{ij}}{r}\right)^6 \right]

The SASA-dependent term modulates the interaction based on the solvent exposure of the residues.

----

* **type**: ``NonBonded``, ``KimHummer``
* **parameters**:

  * ``cutOffNPFactor``: ``real``: Interaction range for LJ as a multiple of sigma
  * ``cutOffDHFactor``: ``real``: Interaction range for DH as a multiple of the Debye length
  * ``dielectricConstant``: ``real``: Relative permittivity of the medium :math:`\epsilon_r`
  * ``debyeLength``: ``real``: Debye length :math:`\lambda_D` :math:`[distance]`
  * ``sasaModel``: ``string``: SASA model to use (options: "A", "B", "C", "D", "E", "F")
  * ``zeroEnergy``: ``real``: Energy shift parameter (default: 0.01) :math:`[energy]`

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
       "zeroEnergy":0.01,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["ALA", "ALA", 1.0, 1.0],
       ["ALA", "GLY", 1.2, 0.9],
       ["GLY", "GLY", 0.8, 1.1]
     ]
   }

.. note::
   The KimHummer potential uses charge, radius, and SASA information stored in the particle data. Ensure that these properties are properly set for each particle before using this potential.

.. warning::
   The KimHummer potential is specifically designed for coarse-grained protein simulations and may not be suitable for other types of systems without modification.
