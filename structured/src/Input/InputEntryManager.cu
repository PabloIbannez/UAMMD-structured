#include "Input/InputEntryManager.cuh"

namespace uammd{
namespace structured{

    void InputEntryManager::loadEntriesInfo(){

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

    InputEntryManager::InputEntryManager(std::shared_ptr<ExtendedSystem> sys,
                                         std::vector<std::string> path):
                                         sys(sys),
                                         path(path){
        this->loadEntriesInfo();
        for(auto entry : entriesInfo){
            System::log<System::DEBUG>("[InputEntryManager] (%s) Entry %s found.",
                                        path.back().c_str(),entry.second.name.c_str());
        }
    }

    bool InputEntryManager::isEntryPresent(std::string name){
        if (entriesInfo.count(name)==0){
            return false;
        }
        return true;
    }

    bool InputEntryManager::isEntrySubTypePresent(std::string entryType,
                                                  std::string entrySubType){
        for(auto e : entriesInfo){
            if(entryType == e.second.entryType and entrySubType == e.second.entrySubType){
                return true;
            }
        }
        return false;
    }

    bool InputEntryManager::isEntryClassPresent(std::string entryType){
        for(auto e : entriesInfo){
            if(entryType == e.second.entryType){
                return true;
            }
        }
        return false;
    }

    InputEntryManager::entryInfo& InputEntryManager::getEntryInfo(std::string name){
        if (!isEntryPresent(name)){
            System::log<System::CRITICAL>("[InputEntryManager] (%s) Requested entry %s not found.",
                                           path.back().c_str(),name.c_str());
        }
        return entriesInfo[name];
    }

    std::map<std::string,InputEntryManager::entryInfo>& InputEntryManager::getEntriesInfo(){
        return entriesInfo;
    }

    std::vector<std::string> InputEntryManager::getEntriesBySubType(std::string entryType,std::string entrySubType){
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

    std::vector<std::string> InputEntryManager::getEntriesByClass(std::string entryType){
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

    void InputEntryManager::checkEntriesUsed(){
        for(auto e : entriesInfo){
            if(!e.second.used){
                System::log<System::WARNING>("[InputEntryManager] (%s) Entry \"%s\" not used.",
                                              path.back().c_str(),e.second.name.c_str());
            }
        }
    }

}}
