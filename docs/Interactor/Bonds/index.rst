Bonds
=====

Bonded interactions are a specific type of interactions that are applied to a predefined set of particles, 
identified by their IDs, and remain constant throughout the simulation. 
These interactions are used to model a wide range of phenomena, 
the most paradigmatic being two particles connected by, for example, 
a harmonic bond. This means that two particles, which we identify using their IDs i and j, 
are interacting with a potential of the type :math:`U(r) = 1/2 K (r - r_0)^2`, where :math:`r = |\vec{r}_{ij}| = |\vec{r}_j - \vec{r}_i|` 
and :math:`K` and :math:`r_0` are the spring constant and the equilibrium distance, 
respectively, which are the characteristic parameters of the potential. 
UAMMD-structured offers a variety of these types of interactions, 
which are categorized based on the number of particles involved in the bond.

----

Available bonds:

.. toctree::
   :maxdepth: 1

   Bond1/index
   Bond2/index
   Bond3/index
   Bond4/index
