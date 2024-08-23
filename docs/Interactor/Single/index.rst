Single
======

In **Single** interaction types, each particle is subjected to a potential. 
This approach allows for the modeling of a variety of scenarios:

.. code-block:: yaml

   topology:
     forceField:
       # ...
       external_example:
         type: ["External", "ConstantForce"]
         parameters:
           constantForce: [0.0, 0.0, 10.0]
       surface_example:
         type: ["Surface", "SurfaceWCAType2"]
         parameters:
           surfacePosition: 0.0
         labels: ["name", "epsilon", "sigma"]
         data:
           - ["A", 1.0, 2.5]
           - ["B", 1.5, 3.0]
       # ...

As demonstrated in the example, such interactions are effective for simulating external fields. 
The given example illustrates the application of a constant force, 
but options also exist for modeling external electric fields or oscillating magnetic fields. 
Additionally, these interactions can incorporate the effects of surfaces. 
UAMMD-structured offers a range of these interaction types.

The implementation process for these interactions is straightforward. In UAMMD-structured,
the implementation is derived from the approach followed in UAMMD, 
tailored to align with input system of UAMMD-structured.

----

Single interactions are categorized into two types:

.. toctree::
   :maxdepth: 1

   External/index
   Surface/index
