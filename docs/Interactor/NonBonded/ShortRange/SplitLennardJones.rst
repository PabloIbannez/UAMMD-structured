SplitLennardJones
-----------------

The SplitLennardJones potential is a modified version of the Lennard-Jones potential that uses different energy scales 
for the repulsive and attractive parts. This allows for more flexibility in tuning the interaction between particles.

.. math::

    U = \begin{cases}
    4\epsilon_r \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right] + (\epsilon_r - \epsilon_a) & \text{if } r < r_c \\
    4\epsilon_a \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right] & \text{if } r_c \leq r < r_{\text{cutoff}}
    \end{cases}

where:

* :math:`\epsilon_r` is the repulsive energy scale
* :math:`\epsilon_a` is the attractive energy scale
* :math:`\sigma` is the distance at which the potential is zero
* :math:`r` is the distance between particles
* :math:`r_c = 2^{1/6}\sigma` is the cutoff distance between repulsive and attractive regimes

----

* **type**: ``NonBonded``, ``SplitLennardJones``
* **parameters**:

  * ``cutOffFactor``: ``real``: Interaction range as a multiple of sigma
  * ``epsilon_r``: ``real``: Repulsive energy scale :math:`[energy]`
  * ``epsilon_a``: ``real``: Attractive energy scale :math:`[energy]`

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Potential well depth :math:`[energy]`
  * ``sigma``: ``real``: Zero-potential distance :math:`[distance]`

Example:

.. code-block::

   "splitLJ":{
     "type":["NonBonded","SplitLennardJones"],
     "parameters":{
       "cutOffFactor":2.5,
       "epsilon_r":1.5,
       "epsilon_a":1.0,
       "condition":"all"
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 1.2, 0.9],
       ["B", "B", 0.8, 1.1]
     ]
   }

.. tip::
   The SplitLennardJones potential allows for independent control of the repulsive and attractive parts of the interaction, which can be useful for modeling systems with complex interparticle interactions.

.. note::
   For an explanation of how select the neighbors list (and the use of the `condition` parameter), see the common documentation for the non-bonded potentials.
