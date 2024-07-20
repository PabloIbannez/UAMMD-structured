#include "GlobalData/Types/TypesHandler.cuh"
#include "GlobalData/Types/TypesFactory.cuh"
#include "GlobalData/Types/Types/Types.cuh"

namespace uammd{
namespace structured{
namespace Types{

    struct Basic_{

        template<typename T>
        static void loadType(std::map<std::string,std::map<std::string,real>>& nameToData,
                             std::map<std::string,T>& typeData){

            std::string name = typeData.at("name");

            nameToData[name]["mass"]   = real(typeData.at("mass"));
            nameToData[name]["radius"] = real(typeData.at("radius"));
            nameToData[name]["charge"] = real(typeData.at("charge"));

            System::log<System::MESSAGE>("[Basic] Loaded type %s, mass: %f, radius: %f, charge: %f",
                                         name.c_str(), nameToData[name]["mass"],
                                         nameToData[name]["radius"], nameToData[name]["charge"]);
        }

        static void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd,
                                              std::map<int,std::string>&    idToName,
                                              std::map<std::string,std::map<std::string,real>>& nameToData){

            int N = pd->getNumParticles();

            auto pos     = pd->getPos(access::location::cpu,access::mode::read);

            //Check if mass,radius or charge are already defined

            bool massDefined   = pd->isMassAllocated();
            bool radiusDefined = pd->isRadiusAllocated();
            bool chargeDefined = pd->isChargeAllocated();

            if(massDefined){
                System::log<System::WARNING>("[Basic] Mass is already defined, ignoring mass from type");
            }
            if(radiusDefined){
                System::log<System::WARNING>("[Basic] Radius is already defined, ignoring radius from type");
            }
            if(chargeDefined){
                System::log<System::WARNING>("[Basic] Charge is already defined, ignoring charge from type");
            }

            auto mass    = pd->getMass(access::location::cpu,   access::mode::write);
            auto radius  = pd->getRadius(access::location::cpu, access::mode::write);
            auto charge  = pd->getCharge(access::location::cpu, access::mode::write);

            for(int i = 0; i < N; i++){

                std::string name = idToName.at(int(pos[i].w));

                if(!massDefined)  { mass[i]   = nameToData[name]["mass"]; }
                if(!radiusDefined){ radius[i] = nameToData[name]["radius"]; }
                if(!chargeDefined){ charge[i] = nameToData[name]["charge"]; }

                System::log<System::DEBUG1>("[Basic] Loading type for particle %d, mass: %f, radius: %f, charge: %f",i,mass[i],radius[i],charge[i]);
            }
        }
    };

    using Basic = Types_<Basic_>;

}}}

REGISTER_TYPES(
    Types,Basic,
    uammd::structured::Types::Basic
)
