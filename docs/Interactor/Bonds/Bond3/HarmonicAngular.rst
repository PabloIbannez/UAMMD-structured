HarmonicAngular
---------------

The HarmonicAngular potential is a simple model for angular bonds that applies a harmonic restraint on the angle between three particles.

.. math::

    U = \frac{1}{2}K(\theta - \theta_0)^2

where:

* :math:`K` is the angular spring constant
* :math:`\theta` is the current angle between the three particles
* :math:`\theta_0` is the equilibrium angle

----

* **type**: ``Bond3``, ``HarmonicAngular``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle
  * ``K``: ``real``: Angular spring constant :math:`[energy]/[angle]^2`
  * ``theta0``: ``real``: Equilibrium angle :math:`[angle]`

Example:

.. code-block::

   "harmonicAngularBonds":{
     "type":["Bond3","HarmonicAngular"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k", "K", "theta0"],
     "data":[[0, 1, 2, 100.0, 1.57],
             [1, 2, 3, 100.0, 1.57],
             [2, 3, 4, 100.0, 1.57]]
   }

HarmonicAngularCommon_K
~~~~~~~~~~~~~~~~~~~~~~~

HarmonicAngular bonds variant with a common angular spring constant (``K``) for all bonds.

----

* **type**: ``Bond3``, ``HarmonicAngularCommon_K``
* **parameters**:

  * ``K``: ``real``: Common angular spring constant for all bonds :math:`[energy]/[angle]^2`

* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle
  * ``theta0``: ``real``: Equilibrium angle :math:`[angle]`

Example:

.. code-block::

   "harmonicAngularBondsCommonK":{
     "type":["Bond3","HarmonicAngularCommon_K"],
     "parameters":{"K":100.0},
     "labels":["id_i", "id_j", "id_k", "theta0"],
     "data":[[0, 1, 2, 1.57],
             [1, 2, 3, 1.57],
             [2, 3, 4, 1.57]]
   }

HarmonicAngularCommon_K_theta0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HarmonicAngular bonds variant with common parameters (``K`` and ``theta0``) for all bonds.

----

* **type**: ``Bond3``, ``HarmonicAngularCommon_K_theta0``
* **parameters**:

  * ``K``: ``real``: Common angular spring constant for all bonds :math:`[energy]/[angle]^2`
  * ``theta0``: ``real``: Common equilibrium angle for all bonds :math:`[angle]`

* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the central particle
  * ``id_k``: ``int``: Id of the third particle

Example:

.. code-block::

   "harmonicAngularBondsCommonKTheta0":{
     "type":["Bond3","HarmonicAngularCommon_K_theta0"],
     "parameters":{"K":100.0,
                   "theta0":1.57},
     "labels":["id_i", "id_j", "id_k"],
     "data":[[0, 1, 2],
             [1, 2, 3],
             [2, 3, 4]]
   }
