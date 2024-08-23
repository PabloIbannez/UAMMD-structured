Gaussian
--------

The Gaussian potential creates a soft, bell-shaped interaction between particles. It's useful for modeling polymer bonds or other soft interactions where a harmonic potential might be too stiff.

.. math::

    U = -E \exp\left(-\frac{(r-r_0)^2}{2D}\right)

where:

* :math:`E` is the well depth
* :math:`r_0` is the equilibrium distance
* :math:`D` is related to the width of the Gaussian

----

* **type**: ``Bond2``, ``Gaussian``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``E``   : ``real``: Well depth :math:`[energy]`
  * ``r0``  : ``real``: Equilibrium distance :math:`[distance]`
  * ``D``   : ``real``: Width parameter :math:`[distance]^2`

Example:

.. code-block::

   "gaussianBonds":{
     "type":["Bond2","Gaussian"],
     "parameters":{},
     "labels":["id_i", "id_j", "E", "r0", "D"],
     "data":[[0, 1, 1.0, 1.0, 0.1],
             [1, 2, 1.0, 1.1, 0.1],
             [2, 3, 1.0, 0.9, 0.1]]
   }

GaussianCommon_E_r0_D
~~~~~~~~~~~~~~~~~~~~~

Gaussian bonds variant with common parameters (``E``, ``r0``, and ``D``) for all bonds.

----

* **type**: ``Bond2``, ``GaussianCommon_E_r0_D``
* **parameters**:

  * ``E`` : ``real``: Common well depth for all bonds :math:`[energy]`
  * ``r0``: ``real``: Common equilibrium distance for all bonds :math:`[distance]`
  * ``D`` : ``real``: Common width parameter for all bonds :math:`[distance]^2`

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle

Example:

.. code-block::

   "gaussianBondsCommon":{
     "type":["Bond2","GaussianCommon_E_r0_D"],
     "parameters":{"E":1.0,
                   "r0":1.0,
                   "D":0.1},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [1, 2],
             [2, 3]]
   }
