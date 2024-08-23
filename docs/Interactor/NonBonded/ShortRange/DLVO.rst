DLVO
----

The DLVO (Derjaguin-Landau-Verwey-Overbeek) potential combines electrostatic repulsion with van der Waals attraction to model colloidal interactions in aqueous solutions. It is the sum of a repulsive electrostatic term (Debye-Hückel) and an attractive van der Waals term (typically modeled as a Lennard-Jones potential).

.. math::

    U_{DLVO} = U_{DH} + U_{LJ}

where:

* :math:`U_{DH}` is the Debye-Hückel potential
* :math:`U_{LJ}` is the Lennard-Jones potential

The Debye-Hückel term is given by:

.. math::

    U_{DH} = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r} \frac{e^{-r/\lambda_D}}{r}

The Lennard-Jones term can be any of the three types (Type1, Type2, or Type3) as described in the LennardJones potential documentation.

----

* **type**: ``NonBonded``, ``DLVOType1`` (or ``DLVOType2`` or ``DLVOType3``)
* **parameters**:

  * ``cutOffNPFactor``: ``real``: Interaction range for LJ as a multiple of sigma
  * ``cutOffDHFactor``: ``real``: Interaction range for DH as a multiple of the Debye length
  * ``dielectricConstant``: ``real``: Relative permittivity of the medium :math:`\epsilon_r`
  * ``debyeLength``: ``real``: Debye length :math:`\lambda_D` :math:`[distance]`

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: LJ potential well depth :math:`[energy]`
  * ``sigma``: ``real``: LJ zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "dlvo":{
     "type":["NonBonded","DLVOType1"],
     "parameters":{
       "cutOffNPFactor":2.5,
       "cutOffDHFactor":3.0,
       "dielectricConstant":78.5,
       "debyeLength":1.0,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 1.2, 0.9],
       ["B", "B", 0.8, 1.1]
     ]
   }

.. note::
   The DLVO potential uses the charge and radius information stored in the particle data. Ensure that these properties are properly set for each particle before using this potential.

.. warning::
   The DLVO theory has limitations and may not accurately describe all colloidal systems, especially at short distances or in highly concentrated solutions.
