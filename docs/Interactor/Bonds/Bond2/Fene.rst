Fene
-----

FENE (Finitely Extensible Nonlinear Elastic) potential to model interactions.
It characterizes a bond that can be stretched to a defined maximum length, beyond which it cannot be extended.

.. math::

    U = -\frac{1}{2}K R_0^2 \ln \left[ 1 - \left( \frac{r-r_0}{R_0} \right)^2 \right]

----

* **type**: ``Bond2``, ``Fene``.
* **parameters**:
  ``None``

* **data**:
  
  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``K``   : ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``r0``  : ``real``: Equilibrium distance :math:`[distance]`
  * ``R0``  : ``real``: Maximum extension of the bond :math:`[distance]`

    
Example:
    
.. code-block::

   "entryName":{
     "type":["Bond2","Fene"],
     "parameters":{},
     "labels":["id_i", "id_j", "r0", "K", "R0"],
     "data":[[0, 1, 1.0, 2.0, 3.0],
             [3, 5, 1.5, 1.0, 2.0]]
   }


FeneCommon_K_R0
~~~~~~~~~~~~~~~

FENE bonds variant with common parameters (``K`` and ``R0``) for all bonds.

----

* **type**: ``Bond2``, ``FeneCommon_K_R0``.
* **parameters**:

  * ``K``  : ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``R0``: ``real``: Maximum extension of the bond :math:`[distance]`

* **data**:
  
  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``r0``: ``real``: Equilibrium distance :math:`[distance]`

----
    
Example:
  
.. code-block::

   "entryName":{
     "type":["Bond2","FeneCommon_K_R0"],
     "parameters":{"K":2.0,
                   "R0:5.0},
     "labels":["id_i", "id_j", "r0"],
     "data":[[0, 1, 3.0],
             [3, 5, 2.0]]
   }

FeneCommon_r0_K_R0
~~~~~~~~~~~~~~~~~~

FENE bonds variant with common parameters (``r0``, ``K`` and ``R0``) for all bonds.

----

* **type**: ``Bond2``, ``FeneCommon_K_R0``.
* **parameters**:

  * ``K``  : ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``R0``: ``real``: Maximum extension of the bond :math:`[distance]`
  * ``r0``: ``real``: Equilibrium distance :math:`[distance]`

* **data**:
  
  * ``id_i``: ``int``: Id of one particle.
  * ``id_j``: ``int``: Id of the other particle.

----

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","FeneCommon_K_R0"],
     "parameters":{"K":2.0,
                   "R0:5.0,
		   "r0":1.0},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [3, 5]]
   }

