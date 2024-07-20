#pragma once

namespace uammd{
namespace structured{
namespace Types{

    template<class T>
    class Types_ : public TypesHandler{

        public:

            Types_(DataEntry& data): TypesHandler(data){
                auto typesData = data.getDataMap();
                for(auto& type: typesData){
                    try{
                        T::loadType(this->nameToData,type);
                    } catch (std::exception& e){
                        System::log<System::CRITICAL>("[Types] Error loading types : %s",e.what());
                    }
                }
            }

            void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd) override {
                T::loadTypesIntoParticleData(pd,this->idToName,this->nameToData);
            }
    };

}}}
