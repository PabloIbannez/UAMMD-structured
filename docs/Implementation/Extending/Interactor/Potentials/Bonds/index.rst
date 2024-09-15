Bonds
=====

Bond Interactors iterates over all bonds the user defines (lines of the data input),
the parameters of the bonds is managed by BondParameters struct, make sure you
include one variable for each label in the input.

Although the name of the structure is ``BondParameters``, the information is read
from the 'data' input. This means that for ``Bonds``, the 'data' input is reserved
for bond-specific information. If you want to parameterize something different
that is not bond-specific information, such as a constant common to all bonds,
you can always do this through the 'parameters' input.

The main implementation difference between Bond1, 2, 3 and 4 is the arguments of computables
functions. AngularBond3 and AngularBond4 are alternatives of Bond3 and Bond4 when the Potential
only deppends on the angle between particles (most of the times).

.. toctree::
   :maxdepth: 4

   Bond1
   Bond2
   Bond3
   Bond4
