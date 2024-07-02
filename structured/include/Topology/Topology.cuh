#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "GlobalData/GlobalData.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "Interactor/InteractorLoader.cuh"

namespace uammd{
namespace structured{

class Topology{

    private:

        std::shared_ptr<ExtendedSystem>      sys;
        std::shared_ptr<GlobalData>           gd;
        std::shared_ptr<ExtendedParticleData> pd;

        std::vector<std::string> path;

        //////////////////////////////

        std::shared_ptr<InputEntryManager>   forceFieldInfo;

        //////////////////////////////

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;
        std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>> VConListSet;
        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> interactors;

        //////////////////////////////

        void loadStructure();

        //Force field

        void loadGroups();
        void loadNeighbourLists();
        void loadInteractors();

    public:

        Topology(std::shared_ptr<ExtendedSystem>       sys,
                 std::shared_ptr<GlobalData>            gd,
                 std::shared_ptr<ExtendedParticleData>  pd,
                 std::vector<std::string> path);

        Topology(std::shared_ptr<ExtendedSystem>       sys,
                 std::shared_ptr<GlobalData>            gd,
                 std::shared_ptr<ExtendedParticleData>  pd);

        //Add a new interactor to the system
        void addInteractor(std::shared_ptr<typename uammd::Interactor> interactor, std::string name);

        //Getters
        std::shared_ptr<ExtendedSystem> getSystem(){ return sys; }

        std::shared_ptr<GlobalData> getGlobalData(){ return gd; }

        std::shared_ptr<ExtendedParticleData> getParticleData(){ return pd; }

        //Get particle group
        std::shared_ptr<ParticleGroup> getParticleGroup(std::string name){
            if(groups.count(name) == 0){
                System::log<System::CRITICAL>("[Topology] Error getting group \"%s\","
                                              " group does not exist.",name.c_str());
            }
            return groups[name];
        }

        //Get default groups, all particles
        std::shared_ptr<ParticleGroup> getParticleGroup(){
            return getParticleGroup("All");
        }

        //Get all groups
        std::map<std::string,std::shared_ptr<ParticleGroup>> getParticleGroups(){
            return groups;
        }

        //Get neighbour list
        std::shared_ptr<VerletConditionalListSetBase> getNeighbourList(std::string name){
            if(VConListSet.count(name) == 0){
                System::log<System::CRITICAL>("[Topology] Error getting neighbour list \"%s\", "
                                              "neighbour list does not exist.",name.c_str());
            }
            return VConListSet[name];
        }

        //Get all neighbour lists
        std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>> getNeighbourLists(){
            return VConListSet;
        }

        //Get interactor
        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> getInteractorsByClass(std::string intClass);

        std::shared_ptr<typename uammd::Interactor> getInteractor(std::string name){
            if(interactors.count(name) == 0){
                System::log<System::CRITICAL>("[Topology] Error getting interactor \"%s\","
                                              " interactor has not been added.",name.c_str());
            }
            return interactors[name];
        }

        //Get all interactors
        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> getInteractors(){
            return interactors;
        }

        InputEntryManager::entryInfo getInteractorInfo(std::string name){
            if(interactors.count(name) == 0){
                System::log<System::CRITICAL>("[Topology] Error getting interactor data \"%s\","
                                              " interactor has not been added.",name.c_str());
            }
            return forceFieldInfo->getEntryInfo(name);
        }

};

}}


