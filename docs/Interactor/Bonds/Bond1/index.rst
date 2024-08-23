Bond1
=====

Often, we want to apply specific potentials to one particle only. 
In these cases, we indicate the ID of the particle. 
This type of bond is encompassed under the type "Bond1":

.. code-block:: yaml

   topology:
     forceField:
       # ...
       bond1_example:
         type: ["Bond1", "FixedHarmonic"]
         labels: ["id_i", "K", "r0", "position"]
         data:
           - [0, 1.5, 0.0, [0.1, 2.3, 5.8]]
           - [1, 2.5, 0.0, [2.2, 0.3, 0.6]]
           - [2, 1.3, 0.0, [1.1, 3.5, 0.9]]
           - [3, 2.1, 0.0, [0.5, 1.1, 1.7]]
       # ...

In the example, certain particles identified by their IDs 0, 1, 2, 3 are subjected to a potential of the type:

.. math::

   U(\mathbf{r}_i) = \frac{1}{2}K(|\mathbf{r}_{i} - \mathbf{pos}| - r_{0})^2

Where **pos** represents a specific position in space. 
This potential is commonly used to fix the position of particles in space. 
In addition to this potential, "Bond1" also allows us to apply forces or torques to specific particles. 
There is also a version of this potential where each term is multiplied by the parameter :math:`\lambda^2` of TI, 
enabling us to use these types of potentials to calculate free energy differences between different particle configurations.

----

Available Bond1 potentials:

.. toctree::
   :maxdepth: 1

   FixedHarmonic
   ConstantForce
   LambdaFixedHarmonicAnisotropic
   FixedHarmonicAnisotropic
