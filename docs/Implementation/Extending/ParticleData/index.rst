Extending ParticleData
======================

UAMMD-structured inherits its ``ExtendedParticleData`` class from
UAMMD's ``ParticleData``, so if you are interested in fully understand how is it
implemented, you should read `UAMMD <https://uammd.readthedocs.io/en/latest/ParticleData.html>`_ documentation.

However, UAMMD-structured offers a very simple way to add a new property
to ParticleData without needing to understand how ParticleData works at
a low level. Simply add the new property to the Data.json file,
and the compiler will automatically define the functions ``pd->getMyProperty()``, ``pd->isMyPropertyAllocated()``, and ``pd->getPropertyIfAllocate()``.
This will allow you to use your new property anywhere in the code.

There is an example of how to add to ``Data.json`` a new property. Write the first string with lowercase and the second with uppercase
is not mandatory, but recommendable to follow the UAMMD standard. ``Data.json`` contains only the properties that have been added
in UAMMD-structured and are not available in UAMMD.

.. code-block::
   :emphasize-lines: 31

    {
        "ParticleData": [
            ["resId", "ResId", "int"],
            ["chainId", "ChainId", "int"],
            ["modelId", "ModelId", "int"],
            ["batchId", "BatchId", "int"],
            ["state", "State", "int4"],
            ["stateInfo", "StateInfo", "real4"],
            ["stress", "Stress", "tensor3"],
            ["frictionConstant", "FrictionConstant", "real"],
            ["TranslationalSelfDiffusion", "TranslationalSelfDiffusion", "real"],
            ["rotationalSelfDiffusion", "RotationalSelfDiffusion", "real"],
            ["SASA", "SASA", "real"],
            ["SASAweight", "SASAweight", "real"],
            ["innerRadius", "InnerRadius", "real"],
            ["epsilon", "Epsilon", "real"],
            ["surface", "Surface", "real"],
            ["patchPos", "PatchPos", "real4"],
            ["patchVector", "PatchVector", "real4"],
            ["parentIndex", "ParentIndex", "int"],
            ["anisotropy", "Anisotropy", "real"],
            ["magneticField", "MagneticField", "real4"],
            ["magnetization", "Magnetization", "real4"],
            ["tentativeState", "TentativeState", "int4"],
            ["transitionProbability", "TransitionProbability", "real"],
            ["selectedId", "SelectedId", "int"],
            ["hessian", "Hessian", "tensor3"],
            ["lambdaDerivative","LambdaDerivative","real"],
            ["polarizability","Polarizability","real"],
            ["pairwiseForce","PairwiseForce","real4"],
            ["myProperty","MyProperty","DataType"]
        ],


Finally, if you want the property to be initialized from the input through `State <../../../State.html>`_,
also add it to the "State" section of the Data.json file.

.. code-block:: json
   :emphasize-lines: 11

    "State": [
        ["velocity", "Vel", "real3"],
        ["direction", "Dir", "real4"],
        ["innerRadius", "InnerRadius", "real"],
        ["magnetization", "Magnetization", "real4"],
        ["anisotropy", "Anisotropy", "real"],
        ["mass", "Mass", "real"],
        ["polarizability","Polarizability","real"],
        ["radius", "Radius", "real"],
        ["charge", "Charge", "real"],
        ["myProperty","MyProperty","DataType"]
    ],

