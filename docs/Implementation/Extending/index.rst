Extending UAMMD-structured
==========================

To add a new module to UAMMD-structured, you must first include it in the **/structured/Components.json** file under the appropriate section. 

.. code-block:: json

    {
      "Class": {
        "Family": [
          ["type", "subtype", "module.cu"],
          ["..."]
        ]
      },
      "..."
    }

Next, create the corresponding file in **/structured/src/Class/Family/Type/module.cu**. 
The content of each **module.cu** file depends on the specific **Class/Family/Type** being implemented. 
However, all modules share a common requirement: to define the module within the module.cu file, 
a **REGISTER** function must be included, depending on the Type being implemented. 
Multiple modules can be implemented in the same **module.cu** file by simply using several **REGISTER** functions.

.. code-block:: cpp

   REGISTER_CLASS(
       Type, SubType,
       uammd::structured::Type::YourModuleObject
   )

----

Available minimal code examples:

.. toctree::
   :maxdepth: 1

   ParticleData/index
   Global/index
   Interactor/index
   Integrator/index
   SimulationStep/index

