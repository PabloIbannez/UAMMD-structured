Generic Interactor
==================

The most simple version of an Interactor:

.. code-block:: cpp

    #pragma once
    #include "uammd.cuh"
    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "Interactor/Interactor.cuh"
    #include "Interactor/InteractorFactory.cuh"
    
    namespace uammd{
    namespace structured{
    namespace Interactor{
    
        class myInteractor: public Interactor{
            public:
                // The constructor must respect this format
                // Add code you want to be executed when the simulation creates the Interactor.
                myInteractor(std::shared_ptr<GlobalData>           gd,
                             std::shared_ptr<ParticleGroup>        pg,
                             DataEntry& data,
                             std::string name): Interactor(pg,"myInteractor: \"" +name+"\""){}

                // Only mandatory function of an Interactor.
                void sum(Computables comp,cudaStream_t st) override {
                // This function is what the Integrator executes every step
                  }
                };
    }}}

    REGISTER_INTERACTOR(
        myType,myInteractor,
        uammd::structured::Interactor::myInteractor
    )


---------

This is an example were we have a basic Interactor that sum a constant force readed from the .json file to the first particle. 

.. code-block:: cpp

    #pragma once
    #include "uammd.cuh"
    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "Interactor/Interactor.cuh"
    #include "Interactor/InteractorFactory.cuh"
    
    namespace uammd{
    namespace structured{
    namespace Interactor{
    
        class myInteractor: public Interactor{
    
            protected:

                std::shared_ptr<GlobalData> gd;

                //you can define your own atributes to use them in sum()
                real3 constant3;

            public:
                myInteractor(std::shared_ptr<GlobalData>           gd,
                             std::shared_ptr<ParticleGroup>        pg,
                             DataEntry& data,
                             std::string name): Interactor(pg,"myInteractor: \"" +name+"\""), gd(gd){

                    // Read input parameters, you will read them from the .json with the name you give them. 
                    // An error will rise if the data from the .json don't match the dataType 
                    constant3 = data.getParameter<real3>("constant3");

                    // Read input data, you will read an entire column from the .json with the label you give them. 
                    // An error will rise if the data from the .json don't match the dataType 
                    std::vector<real> someInput   = data.getData<real>("someInputLabel");
                    std::vector<real3> someInput3 = data.getData<real3>("someInput3Label");
                }
    
                // Only mandatory function of an Interactor is sum, override it from Interactor class.
                void sum(Computables comp,cudaStream_t st) override {

                    std::shared_ptr<ParticleData> pd = pg->getParticleData(); //
                    System::log<System::MESSAGE>("Computing interaction"); //Output logs
                    if(comp.force){
                      // Sum forces to each particle
                      // For instance, adding a force constant3
                      // of the first particle
                      auto forces = pd->getForce(access::cpu, access::write);
                      forces[0] += constant3;
                    }

                    if(comp.energy){
                      //Sum energies to each particle
                    }

                    // Add any computable
                  }
                };
    }}}

    REGISTER_INTERACTOR(
        myType,myInteractor,
        uammd::structured::Interactor::myInteractor
    )


---------
The .json input file that reads this Interactor should be:

.. code-block::

   "name":{
     "type":["myType","myInteractor"],
     "parameters":{"constant3": [1.0,4.0,-1.0]},
     "labels":["someInput", "someInput3"],
     "data":[[0.0, [0.0,1.0,2.0]],
             [1.0, [3.0,4.0,5.0]],
             [2.0, [6.0,7.0,8.0]],
             .....
             [N, [3N,3N+1,3N+2]]
   }

To register your own Interactor create the file
``src/Interactor/Family/myType/myInteractor.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Interactor":
        "Family":[
                 ["..."],
                 ["myType","myInteractor","myInteractor.cu"]
                 ]
   }

