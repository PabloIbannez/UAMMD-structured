Types
=====

If you want to **group** particles based on your own criterion you must implement a new type.
Adding a new type essentially means choosing which ParticleData properties you want to use
to group your particles. The most basic example is if you want your particle selection
to have, in addition to mass, charge, and radius, your own ``myProperty``.

If you are looking to implement a new ``Property``, go to `ParticleData <../ParticleData.html>`_

.. code-block:: cpp
   :emphasize-lines: 20,27,43,54,55,56,61,70,72

    #include "GlobalData/Types/TypesHandler.cuh"
    #include "GlobalData/Types/TypesFactory.cuh"
    #include "GlobalData/Types/Types/Types.cuh"

    namespace uammd{
    namespace structured{
    namespace Types{

        struct MyType_{

            template<typename T>
            static void loadType(std::map<std::string,std::map<std::string,real>>& nameToData,
                                 std::map<std::string,T>& typeData){

                std::string name = typeData.at("name");

                nameToData[name]["mass"]       = real(typeData.at("mass"));
                nameToData[name]["radius"]     = real(typeData.at("radius"));
                nameToData[name]["charge"]     = real(typeData.at("charge"));
                nameToData[name]["myProperty"] = real(typeData.at("myProperty"));

                System::log<System::MESSAGE>("[MyType] Loaded type %s, mass: %f, radius: %f, charge: %f, myProperty: %f",
                                             name.c_str(),
                                             nameToData[name]["mass"],
                                             nameToData[name]["radius"],
                                             nameToData[name]["charge"],
                                             nameToData[name]["myProperty"]);
            }

            static void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd,
                                                  std::map<int,std::string>&    idToName,
                                                  std::map<std::string,std::map<std::string,real>>& nameToData){

                int N = pd->getNumParticles();

                auto pos     = pd->getPos(access::location::cpu,access::mode::read);

                //Check if mass,radius or charge are already defined

                bool massDefined       = pd->isMassAllocated();
                bool radiusDefined     = pd->isRadiusAllocated();
                bool chargeDefined     = pd->isChargeAllocated();
                bool myPropertyDefined = pd->isMyPropertyAllocated();

                if(massDefined){
                    System::log<System::WARNING>("[MyType] Mass is already defined, ignoring mass from type");
                }
                if(radiusDefined){
                    System::log<System::WARNING>("[MyType] Radius is already defined, ignoring radius from type");
                }
                if(chargeDefined){
                    System::log<System::WARNING>("[MyType] Charge is already defined, ignoring charge from type");
                }
                if(myPropertyDefined){
                    System::log<System::WARNING>("[MyType] myProperty is already defined, ignoring myProperty from type");
                }

                auto mass       = pd->getMass(access::location::cpu,   access::mode::write);
                auto radius     = pd->getRadius(access::location::cpu, access::mode::write);
                auto charge     = pd->getCharge(access::location::cpu, access::mode::write);
                auto myProperty = pd->getMyProperty(access::location::cpu, access::mode::write);

                for(int i = 0; i < N; i++){

                    std::string name = idToName.at(int(pos[i].w));

                    if(!massDefined)      { mass[i]       = nameToData[name]["mass"]; }
                    if(!radiusDefined)    { radius[i]     = nameToData[name]["radius"]; }
                    if(!chargeDefined)    { charge[i]     = nameToData[name]["charge"]; }
                    if(!myPropertyDefined){ myProperty[i] = nameToData[name]["myProperty"]; }

                    System::log<System::DEBUG1>("[MyType] Loading type for particle %d, mass: %f, radius: %f, charge: %f, myProperty: %f",i,mass[i],radius[i],charge[i],myProperty[i]);
                }
            }
        };

        using MyType = Types_<MyType_>;

    }}}

    REGISTER_TYPES(
        Types,MyType,
        uammd::structured::Types::MyType
    )

To register your own Types system create the file
``src/GlobalData/Types/Types/myTypes.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "GlobalData":
        "Types":[
            ["..."],
            ["Types","myTypes","myTypes.cu"]
            ]
   }

