#pragma once

#include "System/ExtendedSystem.cuh"

namespace uammd{
namespace structured{

    class InputEntryManager{

        public:

            struct entryInfo{

                std::string name;

                std::vector<std::string> path;

                std::string entryType;
                std::string entrySubType;

                bool used = false;

            };

        private:

            std::shared_ptr<ExtendedSystem>      sys;

            std::vector<std::string> path;

            //////////////////////////////

            std::map<std::string,entryInfo> entriesInfo;

            //////////////////////////////

            void loadEntriesInfo(){

                auto input = sys->
                             getInput();

                if(input->isEntryAdded(this->path)){
                    std::vector<std::string> entriesName = input->getEntriesList(this->path);

                    for(auto& entry : entriesName){
                        std::vector<std::string> entryPath = this->path;
                        entryPath.push_back(entry);

                        auto typesInfo = input->getTypeEntry(entryPath);
                        std::string entryType    = typesInfo.getType();
                        std::string entrySubType = typesInfo.getSubType();

                        entryInfo entryBuffer;

                        entryBuffer.name = entry;
                        entryBuffer.path = entryPath;
                        entryBuffer.entryType     = entryType;
                        entryBuffer.entrySubType = entrySubType;

                        if(entriesInfo.count(entry) == 0){
                            entriesInfo[entry] = entryBuffer;
                        } else {
                            System::log<System::CRITICAL>("[InputEntryManager] (%s) Entry %s already exists.",
                                                           path.back().c_str(),entry.c_str());
                        }
                    }
                } else {
                    System::log<System::WARNING>("[InputEntryManager] (%s) Entry %s not found.",
                                                  path.back().c_str(),path.back().c_str());
                }
            }

        public:

            InputEntryManager(std::shared_ptr<ExtendedSystem> sys,
                              std::vector<std::string> path):
                              sys(sys),
                              path(path){
                this->loadEntriesInfo();
                for(auto entry : entriesInfo){
                    System::log<System::DEBUG>("[InputEntryManager] (%s) Entry %s found.",
                                                path.back().c_str(),entry.second.name.c_str());
                }
            }

            bool isEntryPresent(std::string name){
                if (entriesInfo.count(name)==0){
                    return false;
                }
                return true;
            }

            entryInfo& getEntryInfo(std::string name){
                if (!isEntryPresent(name)){
                    System::log<System::CRITICAL>("[InputEntryManager] (%s) Requested entry %s not found.",
                                                   path.back().c_str(),name.c_str());
                }
                return entriesInfo[name];
            }

            std::map<std::string,entryInfo>& getEntriesInfo(){
                return entriesInfo;
            }

            bool isEntrySubTypePresent(std::string entryType,std::string entrySubType){
                for(auto e : entriesInfo){
                    if(entryType == e.second.entryType and entrySubType == e.second.entrySubType){
                        return true;
                    }
                }
                return false;
            }

            std::vector<std::string> getEntriesBySubType(std::string entryType,std::string entrySubType){
                if (!isEntrySubTypePresent(entryType,entrySubType)){
                    System::log<System::CRITICAL>("[InputEntryManager] (%s) Requested entry class %s and subType %s, but not found.",
                                                   path.back().c_str(),entryType.c_str(),entrySubType.c_str());
                }

                std::vector<std::string> entriesByClass;
                for(auto e : entriesInfo){
                    if(entryType == e.second.entryType and entrySubType == e.second.entrySubType){
                        entriesByClass.push_back(e.second.name);
                    }
                }

                return entriesByClass;
            }

            bool isEntryClassPresent(std::string entryType){
                for(auto e : entriesInfo){
                    if(entryType == e.second.entryType){
                        return true;
                    }
                }
                return false;
            }

            std::vector<std::string> getEntriesByClass(std::string entryType){
                if (!isEntryClassPresent(entryType)){
                    System::log<System::CRITICAL>("[InputEntryManager] (%s) Requested entry class %s , but not found.",
                                                   path.back().c_str(),entryType.c_str());
                }

                std::vector<std::string> entriesByClass;
                for(auto e : entriesInfo){
                    if(entryType == e.second.entryType){
                        entriesByClass.push_back(e.second.name);
                    }
                }

                return entriesByClass;
            }

            void checkEntriesUsed(){
                for(auto e : entriesInfo){
                    if(!e.second.used){
                        System::log<System::WARNING>("[InputEntryManager] (%s) Entry \"%s\" not used.",
                                                      path.back().c_str(),e.second.name.c_str());
                    }
                }
            }

    };

}}
