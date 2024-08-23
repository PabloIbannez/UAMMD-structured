Set
===
(Many Body Bonded)

For many-body bonded interactions, we refer to a class of potentials that depend not on the individual variables of certain particles, **but on collective variables**. A typical example is a harmonic potential that depends on the distance between the centers of mass of two groups of particles:

.. math::

   U(R) = \frac{1}{2}K(R - R_0)^2

Here, :math:`R = |\vec{dR_12}| = |\vec{R}_2 - \vec{R}_1|` is the distance between the centers of mass of particle groups 1 and 2. Such interactions are commonly used in two scenarios. Firstly, for algorithms that compute free energies along a reaction coordinate, like in umbrella sampling for calculating the Potential of Mean Force (PMF) of two structures based on the distance between their centers of mass. This distance would be the reaction coordinate. Here, a potential like the above is applied at multiple points along the reaction coordinate, varying the R_0 parameter. Techniques like WHAM can then be used to derive the PMF as a function of R. Another application of these potentials is in simulating processes like pulling, where a force is exerted on two groups of particles in opposite directions. Additionally, these potentials can be used to apply torque on groups of particles, similar to the stresses applied by optical tweezers on proteins, DNA, etc. Adding such potentials to a UAMMD-structured simulation is done as presented in the following example:

.. code-block:: yaml

   topology:
     forceField:
       mbb1_example:
         type: ["ManyBodyBond1", "FixedHarmonicCenterOfMass"]
         labels: ["idSet_i", "K", "r0", "position"]
         data:
           - [[0,1,2,3,4,5,6,7,8,9,10], 20.0, 0.0, [0.0,0.0,50.0]]
           - ...
       mbb2_example:
         type: ["ManyBodyBond2", "HarmonicBondBetweenCentersOfMass"]
         labels: ["idSet_i", "idSet_j", "K", "r0"]
         data:
           - [[11,12,13,14,15], [16,17,18,19,20], 10.0, 15.0]
           - ...

These potentials fall under two categories: "ManyBodyBond1" and "ManyBodyBond2". The first category includes potentials applied to a single group of particles, while the second involves potentials applied to two groups of particles [Currently, only these two categories exist, but the implementation can be easily generalized for potentials involving three or more groups of particles]. In the first case, we specify that a potential, which fixes the position of the center of mass, is applied to the particles from indices 0 to 10. The second case corresponds to the equation above, where we indicate that a harmonic potential is applied to the distance between the center of mass of the particle group from indices 11 to 15 and the center of mass of the group from indices 16 to 20. Although in the example only one entry is added under "data" for each potential, an indefinite number can be included. There is just one restriction: for the same potential, **a particle can only be in one group**. 

The complete list of potentials and their parameters can be found in:

.. toctree::
   :maxdepth: 1

   Set1/index
   Set2/index
