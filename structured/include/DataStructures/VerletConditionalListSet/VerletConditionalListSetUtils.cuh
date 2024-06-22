#ifndef __VERLET_CONDITIONAL_LIST_SET_UTILS__
#define __VERLET_CONDITIONAL_LIST_SET_UTILS__

namespace uammd{
namespace structured{
namespace VerletConditionalListSetUtils{

    inline
    std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>> loadVerletConditionalListSetsFromInputEntries(std::shared_ptr<ExtendedSystem>        sys,
                                                                                                                      std::shared_ptr<GlobalData>             gd,
                                                                                                                      std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                                                                                                      std::shared_ptr<InputEntryManager> entries){

        std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>> neigLists;

        for(auto& entry : entries->getEntriesInfo()){
            if(entry.second.entryType == "VerletConditionalListSet"){

                std::shared_ptr<uammd::structured::VerletConditionalListSetBase> nl = VerletConditionalListSetLoaders::loadVerletConditionalListSet(sys,gd,groups,entry.second.path);

                //Check if some nl is already present
                if(neigLists.count(entry.second.name) == 0){
                    neigLists[entry.second.name] = nl;
                    entry.second.used = true;
                } else {
                    System::log<System::CRITICAL>("[VerletConditionalListSetsListFromInputEntries] Error loading Verlet conditional list set,"
                                                  "Verlet conditional list set \"%s\" has already been added.",entry.first.c_str());
                }

            }
        }

        //Print information about neig lists
        for(auto nl : neigLists){
            System::log<System::MESSAGE>("[VerletConditionalListSetsListFromInputEntries] Added Verlet conditional list set \"%s\".",
                                         nl.first.c_str());
        }

        return neigLists;
    }

    inline
    std::shared_ptr<VerletConditionalListSetBase> getNeighbourListFromNeighbourListsList(std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                                                                                         DataEntry& data){
        if(nls.size() == 0){
            System::log<System::CRITICAL>("[NeighbourListFromNeighbourListsList] No neighbour lists have been added");
        }

        if(data.isParameterAdded("neighbourList")){
            std::string nlName = data.getParameter<std::string>("neighbourList");
            if(nls.find(nlName) == nls.end()){
                System::log<System::CRITICAL>("[NeighbourListFromNeighbourListsList] Neighbour list '%s' does not exist",nlName.c_str());
            } else {
                System::log<System::MESSAGE>("[NeighbourListFromNeighbourListsList] Using neighbour list '%s'",nlName.c_str());
                return nls[nlName];
            }
        }

        if(nls.size() > 1){
            System::log<System::CRITICAL>("[NeighbourListFromNeighbourListsList] More than one neighbour list found, please specify which one to use");
        }

        return nls.begin()->second;
    }
}}}
#endif
