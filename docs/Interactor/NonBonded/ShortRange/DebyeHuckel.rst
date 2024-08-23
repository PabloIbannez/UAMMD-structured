DebyeHuckel
-----------

The Debye-Hückel potential models electrostatic interactions between charged particles in an electrolyte solution. It's an approximation of the Poisson-Boltzmann equation, accounting for charge screening effects in ionic solutions.

.. math::

    U = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r} \frac{e^{-r/\lambda_D}}{r}

where:

* :math:`q_i, q_j` are the charges of the interacting particles
* :math:`\epsilon_0` is the vacuum permittivity
* :math:`\epsilon_r` is the relative permittivity of the medium
* :math:`r` is the distance between particles
* :math:`\lambda_D` is the Debye length, characterizing the screening effect

----

* **type**: ``NonBonded``, ``DebyeHuckel``
* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of the Debye length
  * ``dielectricConstant``: ``real``: Relative permittivity of the medium :math:`\epsilon_r`
  * ``debyeLength``: ``real``: Debye length :math:`\lambda_D` :math:`[distance]`

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j

Example:

.. code-block::

   "debyeHuckel":{
     "type":["NonBonded","DebyeHuckel"],
     "parameters":{
       "cutOffFactor":3.0,
       "dielectricConstant":78.5,
       "debyeLength":1.0,
       "condition":"all"
     }
   }

.. note::
   The Debye-Hückel potential uses the charge information stored in the particle data. Ensure that the charges are properly set for each particle before using this potential. The potential strength is automatically calculated based on the charges and the specified dielectric constant.

.. warning::
   The Debye-Hückel potential is an approximation valid for weak electrolyte solutions. For strong electrolytes or highly charged particles, more sophisticated models may be necessary.
