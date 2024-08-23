LennardJonesGaussian
--------------------

The LennardJonesGaussian potential combines a WCA (Weeks-Chandler-Andersen) potential with a Gaussian well, providing a modified interaction that can be useful for modeling certain types of bonds or non-bonded interactions.

.. math::

    U = U_{WCA}(r) + U_{Gaussian}(r)

where:

.. math::

    U_{WCA}(r) = \begin{cases}
        \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 2\left(\frac{\sigma}{r}\right)^6 \right] + \epsilon, & \text{if } r \leq \sigma \\
        0, & \text{if } r > \sigma
    \end{cases}

.. math::

    U_{Gaussian}(r) = -\epsilon\exp\left(-\frac{(r-r_0)^2}{2D}\right)

Here:

* :math:`\epsilon` is the depth of the WCA potential well and the scale of the Gaussian well
* :math:`\sigma` is the distance at which the WCA potential is zero
* :math:`r` is the distance between particles
* :math:`r_0` is the position of the Gaussian well minimum
* :math:`D` is related to the width of the Gaussian well

----

* **type**: ``Bond2``, ``LennardJonesGaussian``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``epsilon``: ``real``: Depth of the WCA potential well and scale of Gaussian well :math:`[energy]`
  * ``sigma``: ``real``: Distance at which the WCA potential is zero :math:`[distance]`
  * ``D``: ``real``: Parameter related to the width of the Gaussian well :math:`[distance]^2`

Example:

.. code-block::

   "ljGaussianBonds":{
     "type":["Bond2","LennardJonesGaussian"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma", "D"],
     "data":[[0, 1, 1.0, 1.0, 0.1],
             [1, 2, 1.0, 1.0, 0.1]]
   }

LennardJonesGaussianCommon_epsilon_D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LennardJonesGaussian bonds variant with common parameters (``epsilon`` and ``D``) for all bonds.

----

* **type**: ``Bond2``, ``LennardJonesGaussianCommon_epsilon_D``
* **parameters**:

  * ``epsilon``: ``real``: Common depth of the WCA potential well and scale of Gaussian well for all bonds :math:`[energy]`
  * ``D``: ``real``: Common parameter related to the width of the Gaussian well for all bonds :math:`[distance]^2`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``sigma``: ``real``: Distance at which the WCA potential is zero :math:`[distance]`

Example:

.. code-block::

   "ljGaussianBondsCommonEpsilonD":{
     "type":["Bond2","LennardJonesGaussianCommon_epsilon_D"],
     "parameters":{"epsilon":1.0,
                   "D":0.1},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 1.0],
             [1, 2, 1.0]]
   }
