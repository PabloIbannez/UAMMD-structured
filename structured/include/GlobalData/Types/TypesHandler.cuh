#pragma once

#include <string>
#include <map>

#include "System/ExtendedSystem.cuh"
#include "ParticleData/ExtendedParticleData.cuh"

namespace uammd{
namespace structured{
namespace Types{

    class TypesHandler{

        protected:

            std::map<int,std::string> idToName;
            std::map<std::string,int> nameToId;

            std::map<std::string,std::map<std::string,real>> nameToData;

        public:

            TypesHandler(DataEntry& data){

                std::vector<std::string> labels = data.getLabels();

                //Check if name in labels
                if(std::find(labels.begin(),labels.end(),"name") == labels.end()){
                    System::log<System::CRITICAL>("[TypesHandler] No name label in data");
                }

                int  typeId=0;
                auto typesData = data.getDataMap();
                for(auto& type : typesData){

                    std::string name = type.at("name");

                    if(nameToId.find(name) != nameToId.end()){
                        System::log<System::CRITICAL>("[TypesHandler] Type name %s already exists",name.c_str());
                    } else {
                        nameToId[name]   = typeId;
                        idToName[typeId] = name;
                        typeId++;
                    }

                }
            }

            int getNumberOfTypes(){
                return idToName.size();
            }

            std::vector<std::string> getTypeNameList(){
                std::vector<std::string> nameList;
                for(auto& name : idToName){
                    nameList.push_back(name.second);
                }
                return nameList;
            }

            std::vector<int> getTypeIdList(){
                std::vector<int> idList;
                for(auto& id : idToName){
                    idList.push_back(id.first);
                }
                return idList;
            }

            std::string getTypeName(int id){
                if(idToName.find(id) == idToName.end()){
                    System::log<System::CRITICAL>("[TypesHandler] Type id %d does not exist",id);
                }
                return idToName[id];
            }

            int getTypeId(std::string name){
                if(nameToId.find(name) == nameToId.end()){
                    System::log<System::CRITICAL>("[TypesHandler] Type name %s does not exist",name.c_str());
                }
                return nameToId[name];
            }

            real getTypeData(std::string name, std::string dataName){
                if(nameToData.find(name) == nameToData.end()){
                    System::log<System::CRITICAL>("[TypesHandler] Type name %s does not exist",name.c_str());
                }
                if(nameToData[name].find(dataName) == nameToData[name].end()){
                    System::log<System::CRITICAL>("[TypesHandler] Type data %s does not exist",dataName.c_str());
                }
                return nameToData[name][dataName];
            }

            real getTypeData(int id, std::string dataName){
                return getTypeData(getTypeName(id),dataName);
            }

            virtual void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd) = 0;
    };


}}}
