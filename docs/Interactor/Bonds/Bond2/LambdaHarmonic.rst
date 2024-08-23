LambdaHarmonic
--------------

The LambdaHarmonic potential is a variant of the Harmonic potential that includes a lambda parameter for smooth switching or alchemical transformations. This is particularly useful in free energy calculations or when gradually introducing or removing bonds.

.. math::

    U = \lambda^n \cdot \frac{1}{2}K(r-r_0)^2

where:

* :math:`\lambda` is the switching parameter (0 ≤ λ ≤ 1)
* :math:`n` is the power to which lambda is raised
* :math:`K` is the spring constant
* :math:`r` is the current distance between the bonded particles
* :math:`r_0` is the equilibrium distance

----

* **type**: ``Bond2``, ``LambdaHarmonic``
* **parameters**:

  * ``n``: ``int``: Power of lambda (default: 2)

* **data**:

  * ``id_i``: ``int``: Id of one particle
  * ``id_j``: ``int``: Id of the other particle
  * ``K``   : ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``r0``  : ``real``: Equilibrium distance :math:`[distance]`

Example:

.. code-block::

   "lambdaHarmonicBonds":{
     "type":["Bond2","LambdaHarmonic"],
     "parameters":{"n":2},
     "labels":["id_i", "id_j", "K", "r0"],
     "data":[[0, 1, 100.0, 1.0],
             [1, 2, 100.0, 1.1],
             [2, 3, 100.0, 0.9]]
   }

.. note::
    The :math:`\lambda` parameter is typically controlled globally for the entire simulation and is not specified per bond.

