========
Topology
========

########
Struture
########

##########
ForceField
##########

**************
DataStructures
**************

VerletConditionalListSet
========================

Condition
---------

All
^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``all``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

Example:

.. code-block:: 

   "entryName":{
     "type":["VerletConditionalListSet","all"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
   }

intramolecular-intermolecular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``intra_inter``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","intra_inter"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
   }

nonExcluded
^^^^^^^^^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``nonExcluded``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**:

   .. list-table::
      :widths: 25 25
      :header-rows: 1
      :align: center

      * - id
        - id_list
      * - ``int``
        - [``int``, ``int`` , ...]

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","nonExcluded"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
     "labels":["id","id_list"],
     "data":[
        [0,[1,2,3]],
        [1,[0,2,3]],
        [2,[0,1,3]],
        [3,[0,1,2]],
        ["...","..."]
     ]
   }

nonExcludedIntramolecular-nonExcludedIntermolecular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Data entry description:

* **type**: ``VerletConditionalListSet``, ``nonExclIntra_nonExclInter``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**:

   .. list-table::
      :widths: 25 25
      :header-rows: 1
      :align: center

      * - id
        - id_list
      * - ``int``
        - [``int``, ``int`` , ...]

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","nonExclIntra_nonExclInter"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
     "labels":["id","id_list"],
     "data":[
        [0,[1,2,3]],
        [1,[0,2,3]],
        [2,[0,1,3]],
        [3,[0,1,2]],
        ["...","..."]
     ]
   }

nonExcludedIntramolecular-nonExcludedIntermolecular-nonExcludedCharged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``nonExclIntra_nonExclInter_nonExclCharged``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**:

   .. list-table::
      :widths: 25 25
      :header-rows: 1
      :align: center

      * - id
        - id_list
      * - ``int``
        - [``int``, ``int`` , ...]

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","nonExclIntra_nonExclInter_nonExclCharged"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
     "labels":["id","id_list"],
     "data":[
        [0,[1,2,3]],
        [1,[0,2,3]],
        [2,[0,1,3]],
        [3,[0,1,2]],
        ["...","..."]
     ]
   }

nonExcludedIntramolecular-nonExcludedIntermolecular-nonExcludedIntermolecularCharged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``nonExclIntra_nonExclInter_nonExclInterCharged``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**:

   .. list-table::
      :widths: 25 25
      :header-rows: 1
      :align: center

      * - id
        - id_list
      * - ``int``
        - [``int``, ``int`` , ...]

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","nonExclIntra_nonExclInter_nonExclInterCharged"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     },
     "labels":["id","id_list"],
     "data":[
        [0,[1,2,3]],
        [1,[0,2,3]],
        [2,[0,1,3]],
        [3,[0,1,2]],
        ["...","..."]
     ]
   }

intermolecularDifferentType
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data entry description:

* **type**: ``VerletConditionalListSet``, ``interDifferentType``.
* **parameters**:

  * ``cutOff`` : ``float``, *optional*, default: 0.

  * ``cutOffVerletFactor`` : ``float``, *optional*, default: 1.1 .

* **data**: ``None``.

Example:

.. code-block::

   "entryName":{
     "type":["VerletConditionalListSet","interDifferentType"]
     "parameters":{
       "cutOff":1.0,
       "cutOffVerletFactor":1.2
     }
   }

**********
Interactor
**********

Bonds
=====

Bonds are lists of particles (ranging from 1 to 4 particles) that interact specifically among themselves through a designated potential.

Bond1
-----

Bond1 represents bonds that act solely on individual particles.

FixedHarmonic
^^^^^^^^^^^^^

This is a harmonic potential that links a particle to a fixed point in space. 

.. math::

     x = x_{particle} - x_{fixed}

.. math::

     y = y_{particle} - y_{fixed}

.. math::
     z = z_{particle} - z_{fixed}

.. math::

     U = \frac{1}{2}K_x(x-x_0)^2 + \frac{1}{2}K_y(y-y_0)^2 + \frac{1}{2}K_z(z-z_0)^2


Data entry description:

* **type**: ``Bond1``, ``FixedHarmonic``.
* **parameters**:
  ``None``

* **data**:
  
  * ``id_i``: ``int``: Id of the particle
  
  * ``K``  : ``real3``: Strength of the bond in each direction.
    
  * ``r_0``: ``real3``: Equilibrium distance in each direction.
    
  * ``position``: ``real3``: Position of the fixed point in the space
    
Example:

.. code-block::

   "entryName":{
     "type":["Bond1","FixedHarmonic"],
     "parameters":{},
     "labels":["id_i", "K", "r_0", "position"],
     "data":[[0, [1, 1, 1], [0, 0, 0], [1, 0, 0]],
             [1, [2, 2, 2], [1, 1, 1], [0, 1, 0]]]
   }

ConstantForce
^^^^^^^^^^^^^

This represents a constant force acting on a particle.

.. math::

    U = -\vec{F}\cdot \vec{r}

The sole parameter is a vector that denotes the applied constant force.

Bond2
-----

Bond2 represents interactions between pairs of particles.

Fene
^^^^


This bond utilizes the FENE (Finitely Extensible Nonlinear Elastic) potential to model interactions.
It characterizes a bond that can be stretched to a defined maximum length, beyond which it cannot be extended.

.. math::

    U = -\frac{1}{2}K R_0^2 \ln \left[ 1 - \left( \frac{r-r_0}{R_0} \right)^2 \right]

Depending on the input parameters it has three sub-classes:

Fene

All the information of the bonds is provieded in the data, each bond can have different constants.

Data entry description:

* **type**: ``Bond2``, ``Fene``.
* **parameters**:
  ``None``

* **data**:
  
  * ``id_i``: ``int``: Id of one particle

  * ``id_j``: ``int``: Id of the other particle
  
  * ``K``   : ``real``: Strength of the bond.
    
  * ``r0``  : ``real``: Equilibrium distance
    
  * ``R0``  : ``real``: Maximum extension of the bond

    
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

The constants K and R0 are parameters common to all the bonds.

Data entry description:

* **type**: ``Bond2``, ``FeneCommon_K_R0``.
* **parameters**:

  * ``K``  : ``real``: Strength of the bond.
    
  * ``R0``: ``real``: Maximum extension of the bond

* **data**:
  
  * ``id_i``: ``int``: Id of one particle

  * ``id_j``: ``int``: Id of the other particle
    
  * ``r0``: ``real``: Equilibrium distance

    
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


The constants K, R0 and r0 are parameters common to all the bonds.

Data entry description:

* **type**: ``Bond2``, ``FeneCommon_K_R0``.
* **parameters**:

  * ``K``  : ``real``: Strength of the bond.
    
  * ``R0``: ``real``: Maximum extension of the bond.

  * ``r0``: ``real``: Equilibrium distance.

* **data**:
  
  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.


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

DebyeHuckel
^^^^^^^^^^^

This bond uses the Debye-Hückel interaction, which describes the electrostatic force between charged particles in a solution,
factoring in the screening effect of surrounding ions.


.. math::

    U = \frac{q_1 q_2}{\epsilon r} \exp \left( -\kappa r \right)


Data entry description:

* **type**: ``Bond2``, ``DebyeHuckel``.
* **parameters**: ``none``

* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.
    
  * ``chgProduct``: ``real``. Product of the charges (q1*q2).

  * ``dielectricConstant``: ``real``. Dielectric constant of the medium.

  * ``debyeLength``: ``real``. Debye length.

  * ``cutOff``: ``real``. Maximum distance at which the interaction is computed.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","DebyeHuckel"],
     "parameters":{},
     "labels":["id_i", "id_j", "chgProduct", "dielectricConstant", "debyeLength", "cutOff"],
     "data":[[0, 1, 0.8, 70.0, 10.2, 100.0],
             [2, 4, 0.8, 70.0, 10.2, 100.0]]
   }

    

Harmonic
^^^^^^^^

This bond uses the standard harmonic interaction between two particles.

.. math::
    U = \frac{1}{2} K (r - r_0)^2

Depending on the input parameters it has three sub-classes:

Harmonic

All the information of the bonds is provieded in the data, each bond can have different constants.

Data entry description:

* **type**: ``Bond2``, ``Harmonic``.
* **parameters**: ``none``

* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.
    
  * ``K``: ``real``. Strength of the bond.

  * ``r0``: ``real``. Equilibirum distance.
    
Example:

.. code-block::

   "entryName":{
     "type":["Bond2","Harmonic"],
     "parameters":{},
     "labels":["id_i", "id_j", "K", "r0"],
     "data":[[0, 1, 10.0, 1.0],
             [2, 4, 20.0, 0.5]]
   }


HarmonicCommon_K

The constant K is common to all the bonds.


Data entry description:

* **type**: ``Bond2``, ``Harmonic``.
* **parameters**:

  * ``K``: ``real``. Strength of the bond.

* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.
  
  * ``r0``: ``real``. Equilibirum distance.
    
Example:

.. code-block::

   "entryName":{
     "type":["Bond2","HarmonicCommon_K"],
     "parameters":{"K":15.0},
     "labels":["id_i", "id_j", "r0"],
     "data":[[0, 1, 1.0],
             [2, 4, 0.5]]
   }


HarmonicCommon_K_r0

The constants K and r0 are parameters common to all the bonds.


Data entry description:

* **type**: ``Bond2``, ``HarmonicCommon_K_r0``.
* **parameters**:

  * ``K``: ``real``. Strength of the bond.

  * ``r0``: ``real``. Equilibirum distance.
      
* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.
      
Example:

.. code-block::

   "entryName":{
     "type":["Bond2","HarmonicCommon_K_r0"],
     "parameters":{"K":15.0,
                  "r0":0.75},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [2, 4]]
   }


Steric
^^^^^^

Purely repulsive bonds, of the form,

.. math::
    U = \epsilon \left( \frac{\sigma}{r} \right)^{n}

Steric6

Steric repulsion with n = 6. The values of r and sigma are given as data of each bond.

Data entry description:

* **type**: ``Bond2``, ``Steric6``.
* **parameters**:``none``
      
* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.

  * ``epsilon``: ``real``. Energy of the bond.

  * ``sigma``: ``real``. Particle diameter
      
Example:

.. code-block::

   "entryName":{
     "type":["Bond2","Steric6"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 2.0, 2.0],
             [2, 4, 3.0, 2.0]]
   }

Steric6Common_epsilon_sigma
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The values of epsilon and sigma are parameters common to all the bonds.

Data entry description:

* **type**: ``Bond2``, ``Steric6Common_epsilon_sigma``.

* **parameters**:

  * ``epsilon``: ``real``. Energy of the bond.

  * ``sigma``: ``real``. Particle diameter.

* **data**:

``id_i``: ``int``. Id of one particle.

``id_j``: ``int``. Id of the other particle.

Example:

.. code-block::

   "entryName":{
   "type":["Bond2","Steric6Common_epsilon_sigma"],
   "parameters":{"epsilon":2.0,
                 "sigma":2.0},
   "labels":["id_i", "id_j"],
   "data":[[0, 1],
   [2, 4]]
   }

Steric12
^^^^^^^^
Steric repulsion with n = 12.

Data entry description:


* **type**: ``Bond2``, ``Steric6Common_epsilon_sigma``.

* **parameters**:``none``

* **data**:

  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.

  * ``epsilon``: ``real``. Energy unit.

  * ``sigma``: ``real``. Distance unit.

Example:

.. code-block::

   "entryName":{
   "type":["Bond2","Steric12"],
   "parameters":{},
   "labels":["id_i", "id_j", "epsilon", "sigma"],
   "data":[[0, 1, 2.5, 1.5],
   [2, 4, 3.5, 2.5]]
   }

Steric12Common_epsilon_sigma
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ϵ and σ are parameters common to all the bonds.


Data entry description:

* **type**: ``Bond2``, ``Steric6Common_epsilon_sigma``.

* **parameters**:

  * ``epsilon``: ``real``. Energy of the bond.

  * ``sigma``: ``real``. Particle diameter.

* **data**:

  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.

Example:

.. code-block::

   "entryName":{
   "type":["Bond2","Steric12Common_epsilon_sigma"],
   "parameters":{"epsilon":2.5,
                 "sigma":1.5},
   "labels":["id_i", "id_j"],
   "data":[[0, 1],
           [2, 4]]
   }

Morse
^^^^^

This bond type is characterized by the Morse potential. The Morse potential describes the energy of the bond as a function of the distance between the atoms, factoring in the equilibrium distance \( r_0 \), the well depth \( D_e \), and the width \( \alpha \) of the potential.

.. math::
    U = D_e \left[ 1 - \exp \left( - \alpha (r - r_0) \right) \right]^2


Morse

In this type, all the parameters of the bond are provided individually for each bond.

Data entry description:

* **type**: ``Bond2``, ``Morse``.
* **parameters**: ``none``

* **data**:
  
  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
  
  * ``E``: ``real``: Depth of the potential well.
    
  * ``r0``: ``real``: Equilibrium distance.

  * ``D``: ``real``: Width parameter of the Morse potential (1/alpha).

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","Morse"],
     "parameters":{},
     "labels":["id_i", "id_j", "E", "r0", "D"],
     "data":[[0, 1, 5.0, 5.0, 2.0],
             [2, 4, 4.5, 5.2, 1.8]]
   }


MorseCommon_D

In this type, the parameter \( D \) is common for all the bonds.

Data entry description:

* **type**: ``Bond2``, ``MorseCommon_D``.
* **parameters**:

  * ``D``: ``real``: Width parameter of the Morse potential (1/alpha).

* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
    
  * ``r0``: ``real``: Equilibrium distance.

  * ``E``: ``real``: Depth of the potential well.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","MorseCommon_D"],
     "parameters":{"D": 1.5},
     "labels":["id_i", "id_j", "r0", "E"],
     "data":[[0, 1, 3.0, 2.0],
             [2, 4, 3.2, 1.8]]
   }


MorseCommon_r0_E_D

In this type, the parameters \( r_0 \), \( D \) and  \( E \) are common for all the bonds.

Data entry description:

* **type**: ``Bond2``, ``MorseCommon_r0_E_D``.
* **parameters**:

  * ``E``: ``real``: Depth of the potential well.
    
  * ``r0``: ``real``: Equilibrium distance.

  * ``D``: ``real``: Width parameter of the Morse potential.
    
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","MorseCommon_r0_E_D"],
     "parameters":{"D": 1.5,
                   "r0": 3.0,
		   "E":2.0},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [2, 4]]
   }

   
MorseWCACommon_eps0

LenardJones
^^^^^^^^^^^

Interactions based on generalizing the Lennard-Jones interactions, with the general form:

.. math::
    U = A \left( \frac{\sigma}{r} \right)^{m} - B \left( \frac{\sigma}{r} \right)^{n}

LennardJonesType1

In this case, A = B = 4, m = 12, and n = 6.


.. math::
    U = 4 \left( \left( \frac{\sigma}{r} \right)^{12} -  \left( \frac{\sigma}{r} \right)^{6} \right)


The values of epsilon and sigma are provided
as data and may differ for each bond.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesType1``.
* **parameters**: ``none``
    
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

  * ``epsilon``: ``real``: Depth of the potential well.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType1"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.5, 2.0],
             [2, 4, 3.0, 1.25]]
   }

   
LennardJonesType1Common_epsilon

Here, A = B = 4, m = 12, and n = 6, but the value of epsilon is provided as a parameter, ensuring
it remains consistent across all bonds.


Data entry description:

* **type**: ``Bond2``, ``LennardJonesType1Common_epsilon``.
* **parameters**:

  * ``epsilon``: ``real``: Depth of the potential well.
  
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
    
  * ``sigma``: ``real``: Size of the particles.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType1Common_epsilon"],
     "parameters":{"epsilon": 2.0},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 2.0],
             [2, 4, 1.25]]
   }


LennardJonesType2

Here A = 1 B = 2, m = 12 and n = 6.



.. math::
    U = \left( \frac{\sigma}{r} \right)^{12} - 2 \left( \frac{\sigma}{r} \right)^{6}

The values of epsilon and sigma are given as data so the might be different in every bond.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesType1``.
* **parameters**: ``none``
    
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

  * ``epsilon``: ``real``: Depth of the potential well.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType1"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.5, 2.0],
             [2, 4, 3.0, 1.25]]
   }

   

LennardJonesType2Common_epsilon


For this type, the energy function is the same as `LennardJonesType2`, but the value of \(ε\)
is provided as a parameter, ensuring it remains consistent across all bonds.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesType2Common_epsilon``.
* **parameters**:

  * ``epsilon``: ``real``: Depth of the potential well.
  
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType2Common_epsilon"],
     "parameters":{"epsilon": 2.0},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 2.0],
             [2, 3, 1.75]]
   }


LennardJonesType3
^^^^^^^^^^^^^^^^^
Here A = 5, B = 6, m = 12 and n = 10.


.. math::
    U = 5 \left( \frac{\sigma}{r} \right)^{12} - 6 \left( \frac{\sigma}{r} \right)^{10}


The values of \(ε\) and \(σ\) are provided as data and might differ for each bond.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesType3``.
* **parameters**: ``none``
    
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

  * ``epsilon``: ``real``: Depth of the potential well.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType3"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.7, 2.2],
             [2, 3, 2.5, 1.5]]
   }

LennardJonesType3Common_epsilon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For this type, the energy function is the same as `LennardJonesType3`, but the value of \(ε\)
is provided as a parameter, ensuring it remains consistent across all bonds.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesType3Common_epsilon``.
* **parameters**:

  * ``epsilon``: ``real``: Depth of the potential well.
  
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesType3Common_epsilon"],
     "parameters":{"epsilon": 1.9},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 2.2],
             [2, 3, 1.5]]
   }


LennardJonesKaranicolasBrooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The energy function for this type includes an additional term that represents the energy required to displace the water molecules in proximity of one particle.

.. math::
    U(r) = \epsilon \left( 13 \left( \frac{\sigma}{r} \right)^{12} - 18 \left( \frac{\sigma}{r} \right)^{10} + 4 \left( \frac{\sigma}{r} \right)^{6} \right)

The values of \(ε\) and \(σ\) are provided as data and might differ for each bond.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesKaranicolasBrooks``.
* **parameters**: ``none``
    
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

  * ``epsilon``: ``real``: Depth of the potential well.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesKaranicolasBrooks"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "sigma"],
     "data":[[0, 1, 1.8, 2.1],
             [2, 3, 2.6, 1.4]]
   }

LennardJonesKaranicolasBrooksCommon_epsilon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For this type, the energy function is the same as `LennardJonesKaranicolasBrooks`, but the value of \(ε\) is provided as a parameter, ensuring it remains consistent across all bonds.

Data entry description:

* **type**: ``Bond2``, ``LennardJonesKaranicolasBrooksCommon_epsilon``.
* **parameters**:

  * ``epsilon``: ``real``: Depth of the potential well.
  
* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
    
  * ``sigma``: ``real``: Size of the particle.


Example:

.. code-block::

   "entryName":{
     "type":["Bond2","LennardJonesKaranicolasBrooksCommon_epsilon"],
     "parameters":{"epsilon": 2.1},
     "labels":["id_i", "id_j", "sigma"],
     "data":[[0, 1, 2.1],
             [2, 3, 1.4]]
   }

LennardJonesGaussian

LennardJonesGaussianCommon_epsilon_D

Gaussian
^^^^^^^^

This bond type is characterized by the Gaussian potential. The Gaussian potential describes the energy of the bond as a function of the distance between the atoms, with a peak at the equilibrium distance \( r_0 \).

.. math::
    U = \epsilon \exp \left( - \frac{(r-r_0)^2}{\sigma^2} \right)


Gaussian
^^^^^^^^

In this type, all the parameters of the bond are provided individually for each bond.

Data entry description:

* **type**: ``Bond2``, ``Gaussian``.
* **parameters**: ``none``

* **data**:
  
  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.
  
  * ``epsilon``: ``real``: Energy depth of the potential.
    
  * ``r0``: ``real``: Equilibrium distance.
  
  * ``sigma``: ``real``: Width of the Gaussian potential.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","Gaussian"],
     "parameters":{},
     "labels":["id_i", "id_j", "epsilon", "r0", "sigma"],
     "data":[[0, 1, 5.0, 1.0, 0.1],
             [2, 4, 4.0, 1.2, 0.15]]
   }


GaussianCommon_E_r0_D
^^^^^^^^^^^^^^^^^^^^^

In this type, the parameters \( \epsilon \), \( r_0 \), and \( \sigma \) are common for all the bonds.

Data entry description:

* **type**: ``Bond2``, ``GaussianCommon_E_r0_D``.
* **parameters**:

  * ``epsilon``: ``real``: Energy depth of the potential.
    
  * ``r0``: ``real``: Equilibrium distance.
  
  * ``sigma``: ``real``: Width of the Gaussian potential.

* **data**:

  * ``id_i``: ``int``: Id of one particle.

  * ``id_j``: ``int``: Id of the other particle.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","GaussianCommon_E_r0_D"],
     "parameters":{"epsilon": 5.0,
                   "r0": 1.0,
                   "sigma": 0.1},
     "labels":["id_i", "id_j"],
     "data":[[0, 1],
             [2, 4]]
   }
   
MaxDistanceRestraint
^^^^^^^^^^^^^^^^^^^^
If the distance between two particles is greater than the maximum distance, the energy is increased by a harmonic potential.

.. math::

   U = \begin{cases}
        0 &\text{if  } r < r_{max} \\
        \frac{1}{2} K (r - r_{max})^2 &\text{if  } r \geq r_{max}
   \end{cases}


Data entry description:

* **type**: ``Bond2``, ``MaxDistanceRestraint``.
* **parameters**: ``none``

* **data**:
  
  * ``id_i``: ``int``. Id of one particle.

  * ``id_j``: ``int``. Id of the other particle.
    
  * ``maxDistance``: ``real``. Distance at which the harmonic bond starts to act.

  * ``K``: ``real``. Strength of the harmonic bond.

Example:

.. code-block::

   "entryName":{
     "type":["Bond2","MaxDistanceRestraing"],
     "parameters":{},
     "labels":["id_i", "id_j", "maxDistance", "K"],
     "data":[[0, 1, 0.8, 30.0, 10.2],
             [2, 4, 0.8, 35.0, 15.5]]
   }


Bond3
-----

Bond3 describes interactions involving triplets of particles.

HarmonicAngular
^^^^^^^^^^^^^^^

The energy of the interaction varies quadratically with the angle formed by the three particles.

.. math::
    U = \frac{1}{2} K (\theta - \theta_0)^2


HarmonicAngular.

The constants \( K \) and \( \theta_0 \) are provided as data entries, allowing them to differ for each bond.

Data entry description:

- **type**: ``Bond3``, ``HarmonicAngular``.
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the central particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``K``: ``real``: Spring constant.
  * ``theta0``: ``real``: Equilibrium angle.

Example:

.. code-block:: json

   {
     "type": ["Bond3","HarmonicAngular"],
     "labels": ["id_i", "id_j", "id_k", "K", "theta_0"],
     "data": [[0, 1, 2, 100.0, 1.570795],
              [2, 3, 4, 95.0,  3.141593]]
   }
   
HarmonicAngularCommon_K

For this type of interaction,  the spring constant \(K\) is provided as a common parameter across all bonds.

Data entry description:

* **type**: ``Bond3``, ``HarmonicAngularCommon_K``.
* **parameters**:
  
  * ``K``: ``real``: Spring constant.

* **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the central particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``theta0``: ``real``: Equilibrium angle.

Example:

.. code-block:: json

   {
     "type": ["Bond3","HarmonicAngularCommon_K"],
     "parameters": {"K": 150.0},
     "labels": ["id_i", "id_j", "id_k", "theta0"],
     "data": [[0, 1, 2, 1.570795]
            , [1, 2, 3, 0.785397]]
   }

HarmonicAngularCommon_K_theta0

For this interaction type, both the spring constant \(K\) and the equilibrium angle \(\theta_0\) are provided as common parameters across all bonds.

Data entry description:

- **type**: ``Bond3``, ``HarmonicAngularCommon_K_theta0``.
- **parameters**:

  * ``K``: ``real``: Spring constant.
  * ``theta_0``: ``real``: Equilibrium angle in degrees.

- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the central particle.
  * ``id_k``: ``int``: Id of the third particle.

Example:

.. code-block:: json

   {
     "type": ["Bond3","HarmonicAngularCommon_K_theta0"],
     "parameters": {"K": 150.0, "theta0": 1.570795},
     "labels": ["id_i", "id_j", "id_k"],
     "data": [[0, 1, 2],
              [1, 2, 3]]
   }

KratkyPorod
^^^^^^^^^^^

This interaction type describes a Kratky-Porod potential based on the cosine of the angle formed by three particles.

.. math::
    U = K \left( 1.0 + \cos \theta \right)

KratkyPorod
  
Data entry description:

- **type**: ``Bond3``, ``KratkyPorod``.
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the central particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``K``: ``real``: Spring constant.

Example:

.. code-block:: json

   {
     "type": ["Bond3","KratkyPorod"],
     "labels": ["id_i", "id_j", "id_k", "K"],
     "data": [[0, 1, 2, 100.0], [2, 3, 4, 95.0]]
   }

KratkyPorodCommon_K
^^^^^^^^^^^^^^^^^^^

For this type of interaction, the spring constant \(K\) is provided as a common parameter across all bonds.

.. math::
    U = K \left( 1.0 + \cos \theta \right)

Data entry description:

- **type**: ``Bond3``, ``KratkyPorodCommon_K``.
- **parameters**:

  * ``K``: ``real``: Common spring.
  
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the central particle.
  * ``id_k``: ``int``: Id of the third particle.

Example:

.. code-block:: json

   {
     "type": ["Bond3","KratkyPorodCommon_K"],
     "parameters": {"K": 100.0},
     "labels": ["id_i", "id_j", "id_k"],
     "data": [[0, 1, 2], [2, 3, 4]]
   }

BestChenHummerAngular
^^^^^^^^^^^^^^^^^^^^^

This interaction type describes the Best-Chen-Hummer angular potential.
It is a double-well potential that can capture the behavior of flexible angles found in certain
molecular systems. The energy function involves two Gaussian distributions with different centers
and amplitudes.


.. math::
   U = -\frac{1}{\gamma} \ln \left[ \exp (-\gamma (K_{\alpha}(\theta - \theta_{\alpha})^2 + \epsilon_{\alpha})) + \exp ( -\gamma K_{\beta}(\theta - \theta_{\beta})^2) \right]


Bond4
-----

Bond4 describes interactions based on the torsion angle formed by four consecutive particles.

Dihedral
^^^^^^^^

Here, the energy of the interaction varies based on the cosine of the torsion angle.

.. math::
    U = K \left( 1.0 + \cos (n \phi + \phi_0) \right)

Dihedral

Data entry description:

- **type**: ``Bond4``, ``Dihedral``.
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the second particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``id_l``: ``int``: Id of the fourth particle.
  * ``K``: ``real``: Spring constant.
  * ``n``: ``int``: Harmonic number.
  * ``phi0``: ``real``: Phase offset.

Example:

.. code-block:: json

   {
     "type": ["Bond4","Dihedral"],
     "parameters":{},
     "labels": ["id_i", "id_j", "id_k", "id_l", "K", "n", "phi0"],
     "data": [[0, 1, 2, 3, 1.0, 2, 3.14159],
              [2, 3, 4, 5, 0.5, 3, 0.0]]
   }
  
DihedralCommon_n_K_phi0

For this type of interaction, the harmonic number \(n\), spring constant \(K\), and phase offset \(\phi_0\) are provided as common parameters across all bonds.

Data entry description:

- **type**: ``Bond4``, ``DihedralCommon_n_K_phi0``.
- **parameters**:

  * ``K``: ``real``: Common spring constant for all the bonds.
  * ``n``: ``int``: Common harmonic number for all the bonds.
  * ``phi0``: ``real``: Common phase offset.
  
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the second particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``id_l``: ``int``: Id of the fourth particle.

Example:

.. code-block:: json

   {
     "type": ["Bond4","DihedralCommon_n_K_phi0"],
     "parameters": {"K": 1.0, "n": 2, "phi_0": 0.0},
     "labels": ["id_i", "id_j", "id_k", "id_l"],
     "data": [[0, 1, 2, 3], [2, 3, 4, 5]]
   }

Dihedral4
^^^^^^^^^

This interaction type captures the energy of dihedral interactions with multiple harmonic terms.

.. math::
    U = \sum_{i=4}^{n} K_i \left( 1.0 + \cos (n_i \phi + \phi_{0_{i}}) \right)

Data entry description:

- **type**: ``Bond4``, ``Dihedral4``.
- **data**:

  * ``id_i``: ``int``: Id of the first particle.
  * ``id_j``: ``int``: Id of the second particle.
  * ``id_k``: ``int``: Id of the third particle.
  * ``id_l``: ``int``: Id of the fourth particle.
  * ``K``: ``array``: List of spring constants.
  * ``phi0``: ``array``: List of phase offsets.

Example:

.. code-block:: json

   {
     "type": ["Bond4","Dihedral4"],
     "labels": ["id_i", "id_j", "id_k", "id_l", "K", "phi0"],
     "data": [[0, 1, 2, 3, [0.7,0.4,0.5,0.56], [0.5,0.4,0.123,0.56]], 
              [2, 3, 4, 5, [0.7,0.4,0.5,0.56], [0.5,0.4,0.123,0.56]]
    }

Single
======

External
--------

ConstantForce
^^^^^^^^^^^^^

.. math::
    U = -F\cdot \vec{r}

SphericalShell
^^^^^^^^^^^^^^

.. math::
    r = \sqrt{(\vec{r}_i - \vec{r}_{shell})^2} - R_{shell}
.. math::
   \sigma = \sigma_{shell} + \sigma_{particle}
.. math::
   U = \epsilon \left( \frac{\sigma}{r} \right)^{12}


ACMagneticField
^^^^^^^^^^^^^^^

ConstantMagneticField
^^^^^^^^^^^^^^^^^^^^^

Surface
-------

LennardJones
^^^^^^^^^^^^

.. math::
    U = A \left( \frac{\sigma}{z} \right)^{m} - B \left( \frac{\sigma}{z} \right)^{n}

SurfaceLennardJonesType1

SurfaceLennardJonesType2

SurfaceWCAType1

SurfaceWCAType2

SurfaceGeneralLennardJonesType1

SurfaceGeneralLennardJonesType2

SurfaceAnchorage
^^^^^^^^^^^^^^^^

.. math::
   d = \frac{z_{surface} - z}{\sigma}

.. math::
   U = \begin{cases}
        0 &\text{if  } z < -\sigma + z_{surface} \\
        \frac{1}{2} \epsilon (d^2 - 1) &\text{if  } -\sigma + z_{surface} \leq z \leq \sigma + z_{surface} \\
        0 &\text{if  } z > \sigma + z_{surface}
   \end{cases}

ParabolaSurface
^^^^^^^^^^^^^^^

Pair
====

NonBonded
---------

DebyeHuckel
^^^^^^^^^^^

Atzberger
^^^^^^^^^

DLVO
^^^^

DLVOType1

DLVOType2

DLVOType3

Clashed
^^^^^^^

KimHummer
^^^^^^^^^

LennardJones
^^^^^^^^^^^^

LennardJonesType1

LennardJonesType2

LennardJonesType3

WCAType1

WCAType2

WCAType3

GeneralLennardJonesType1

GeneralLennardJonesType2

GeneralLennardJonesType3

Steric
^^^^^^

Steric6

Steric12

StericInner6

StericInner12

DipolarMagnetic
^^^^^^^^^^^^^^^


Set
===

Set1
----

ConstantForceOverCenterOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ConstantTorqueOverCenterOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FixedHarmonicCenterOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^

FixedHarmonicCommon_K_r0CenterOfMass

Set2
----

ConstantForceBetweenCentersOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ConstantTorqueBetweenCentersOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

HarmonicBondBetweenCentersOfMass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PatchyParticles
===============

PatchyParticles
---------------

DynamicallyBondedPatchyParticles
--------------------------------

Patches
=======

SurfacePatches
--------------

Linker
^^^^^^

NonBondedPatches
----------------

DistanceSwitch
^^^^^^^^^^^^^^

Helix
^^^^^

Helix2States
^^^^^^^^^^^^

AFM
===

SphericalTip
------------

SphericallyBluntedConicTip
--------------------------

