FixedHarmonic
-------------

:file:`Interactor/Bonds/Bond1/FixedHarmonic.cu`

The FixedHarmonic potential applies a harmonic restraint between a particle and a fixed point in space.

.. math::

    U = \frac{1}{2}K(r - r_0)^2

where:

* :math:`K` is the spring constant
* :math:`r` is the distance between the particle and the fixed point
* :math:`r_0` is the equilibrium distance

----

* **type**: ``Bond1``, ``FixedHarmonic``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the particle
  * ``K``: ``real``: Spring constant :math:`[energy]/[distance]^2`
  * ``r0``: ``real``: Equilibrium distance :math:`[distance]`
  * ``position``: ``real3``: Fixed point coordinates :math:`[distance]`

Example:

.. code-block::

   "fixedHarmonicBonds":{
     "type":["Bond1","FixedHarmonic"],
     "parameters":{},
     "labels":["id_i", "K", "r0", "position"],
     "data":[[0, 100.0, 1.0, [0.0, 0.0, 0.0]],
             [1, 100.0, 1.0, [1.0, 1.0, 1.0]]]
   }
