Bond2
=====

Possibly the most common bonds within this category are those involving two particles. 
These types of interactions are frequently used in modeling systems:

.. code-block:: yaml

   topology:
     forceField:
       # ...
       bond2_example1:
         type: ["Bond2", "Harmonic"]
         labels: ["id_i", "id_j", "K", "r0"]
         data:
           - [5, 6, 2.1, 0.2]
           - [6, 7, 2.6, 1.2]
           - [7, 8, 3.2, 3.8]
       bond2_example2:
         type: ["Bond2", "FeneCommon_K_R0"]
         parameters:
           K: 3.2
           R0: 4.0
         labels: ["id_i", "id_j", "r0"]
         data:
           - [8, 9, 11.3]
           - [9, 10, 1.5]
           - [10, 11, 2.8]
       # ...

The example demonstrates the application of both a harmonic potential and a FENE potential. 
In the latter case, the "FeneCommon_K_R0" version uses shared parameters K and :math:`R_0` for all bonds, 
specified in the "parameters" field. Such variations of potentials are common in UAMMD-structured, 
aiming to simplify the input.

To name a few, Morse potentials are available, as well as Lennard-Jones potentials used in coarse-grained models for proteins. 
This context also includes potentials for mimicking the effects of hydrophobic interactions. 
Moreover, there are versions of the harmonic potential modified with the :math:`\lambda` parameter, enabling their use in Thermodynamic Integration.

----

Available Bond2 potentials:

.. toctree::
   :maxdepth: 1

   Harmonic
   HarmonicAnisotropic
   LambdaHarmonic
   Fene
   Morse
   Gaussian
   Steric
   LennardJones
   LennardJonesKaranicolasBrooks
   LennardJonesGaussian
   LennardJonesSoftCore
   DebyeHuckel
   MaxDistanceRestraint
   HelixExponential
   HelixCosine
