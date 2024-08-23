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

* **type**: ``Bond2``, ``DebyeHuckel``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``chgProduct``: ``real``: Product of charges :math:`q_i q_j`
  * ``dielectricConstant``: ``real``: Dielectric constant of the medium :math:`\epsilon_r`
  * ``debyeLength``: ``real``: Debye length :math:`\lambda_D` :math:`[distance]`
  * ``cutOff``: ``real``: Cutoff distance for the interaction :math:`[distance]`

Example:

.. code-block::

   "debyeHuckelBonds":{
     "type":["Bond2","DebyeHuckel"],
     "parameters":{},
     "labels":["id_i", "id_j", "chgProduct", "dielectricConstant", "debyeLength", "cutOff"],
     "data":[[0, 1, -1.0, 78.5, 1.0, 10.0],
             [1, 2,  1.0, 78.5, 1.0, 10.0],
             [2, 3, -1.0, 78.5, 1.0, 10.0]]
   }

.. warning::
    In this potential, the charge product value is read from the potential data, and the charges of the particles are not used.
