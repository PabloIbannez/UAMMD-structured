BestChenHummerAngular
---------------------

The BestChenHummerAngular potential is a specialized angular potential developed by Best, Chen, and Hummer for protein modeling [BestChenHummer]_. 
It combines two Gaussian-like terms to create a potential with two minima.

.. math::

    U = -\frac{1}{\gamma} \ln \left( e^{-\gamma(k_\alpha(\theta-\theta_\alpha)^2 + \epsilon_\alpha)} + e^{-\gamma k_\beta(\theta-\theta_\beta)^2} \right)

where:

* :math:`\gamma`, :math:`k_\alpha`, :math:`k_\beta`, :math:`\theta_\alpha`, :math:`\theta_\beta`, and :math:`\epsilon_\alpha` are constants defined in the potential
* :math:`\theta` is the current angle between the three particles

----

* **type**: ``Bond3``, ``BestChenHummerAngular``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle

Example:

.. code-block::

   "bestChenHummerAngularBonds":{
     "type":["Bond3","BestChenHummerAngular"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k"],
     "data":[[0, 1, 2],
             [1, 2, 3],
             [2, 3, 4]]
   }

.. warning::

    The potential parameters are hardcoded in the potential implementation (taken from [BestChenHummer]_). 
    The potential is only defined for KcalMol_A units, if other units are selected an error will be raised.

.. [BestChenHummer] Best, Robert B., Yng-Gwei Chen, and Gerhard Hummer. "Slow protein conformational dynamics from multiple experimental structures: the helix/sheet transition of arc repressor." Structure 13, no. 12 (2005): 1755-1763.
