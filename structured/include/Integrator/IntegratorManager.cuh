#pragma once

#include"Integrator/IntegratorBase.cuh"
#include"Integrator/IntegratorLoaders.cuh"

namespace uammd{
namespace structured{

class IntegratorManager{

    private:

        std::shared_ptr<ExtendedSystem>       sys;
        std::shared_ptr<GlobalData>            gd;
        std::shared_ptr<ExtendedParticleData>  pd;

        std::vector<std::string> path;

        //////////////////////////////////////////

        std::shared_ptr<InputEntryManager>   integratorsInfo;

        //////////////////////////////////////////

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;
        std::map<std::string,std::shared_ptr<uammd::Integrator>> integrators;

        struct stepsInfo{
            std::string name;
            uint order;
            ullint steps;
        };

        std::map<std::string,stepsInfo> integratorSteps;

        //////////////////////////////////////////

        void loadGroups(){
            groups = GroupUtils::loadGroupsListFromInputEntries(sys,gd,pd,integratorsInfo);
        }

        void loadSchedule(){

            bool found = false;
            for(auto& entry : integratorsInfo->getEntriesInfo()){
                if(entry.second.entryType == "Schedule" and entry.second.entrySubType == "Integrator"){
                    //Check if found more than one schedule
                    if(found){
                        System::log<System::CRITICAL>("[IntegratorManager] (%s) Multiple schedules found."
                                                      " Only one schedule is allowed per simulation.",path.back().c_str());
                    }
                    found = true;
                    entry.second.used = true;

                    //Load schedule

                    DataEntry data = sys->getInput()->getDataEntry(entry.second.path);
                    for(auto& d : data.getDataMap()){

                        std::string integratorName = d.at("integrator");

                        stepsInfo info;
                        info.name  = integratorName;
                        info.order = d.at("order");
                        //Check if the steps key is present
                        if(d.find("steps") != d.end()){
                            info.steps = d.at("steps");
                        }else{
                            info.steps = std::numeric_limits<ullint>::max();
                        }

                        integratorSteps[integratorName] = info;

                        if(info.steps != std::numeric_limits<ullint>::max()){
                            System::log<System::MESSAGE>("[IntegratorManager] (%s) Integrator '%s' will run for %llu steps with order %u.",
                                                         path.back().c_str(),integratorName.c_str(),info.steps,info.order);
                        } else {
                            System::log<System::WARNING>("[IntegratorManager] (%s) Integrator '%s' will run for an unlimited number of steps with order %u.",
                                                         path.back().c_str(),integratorName.c_str(),info.order);
                        }

                    }
                }
            }

            if(not found){
                System::log<System::CRITICAL>("[IntegratorManager] (%s) No schedule found."
                                              " At least one schedule is required per simulation.",path.back().c_str());
            }

        }

        void loadIntegrators(){

            for(auto& entry : integratorsInfo->getEntriesInfo()){
                if(IntegratorLoader::isIntegratorAvailable(sys,entry.second.path)){
                    std::shared_ptr<uammd::Integrator> integrator = IntegratorLoader::loadIntegrator(sys,gd,groups,entry.second.path);

                    if(integrators.count(entry.second.name) == 0){
                        //Check if integrator is in schedule
                        if(integratorSteps.count(entry.second.name) == 0){
                            System::log<System::CRITICAL>("[IntegratorManager] (%s) Integrator '%s' is not in schedule."
                                                          " All integrators must be in schedule.",path.back().c_str(),entry.second.name.c_str());
                        }

                        integrators[entry.second.name] = integrator;
                        entry.second.used = true;
                    }
                    else{
                        System::log<System::CRITICAL>("[IntegratorManager] (%s) Error loading integrators,"
                                                      "integrator \"%s\" has already been loaded.",path.back().c_str(),entry.second.name.c_str());
                    }

                }
            }

            //Print information about loaded integrators
            for(auto& integrator : integrators){
                System::log<System::MESSAGE>("[IntegratorManager] (%s) Integrator '%s' loaded.",
                                             path.back().c_str(),integrator.first.c_str());
            }
        }

    public:

        IntegratorManager(std::shared_ptr<ExtendedSystem>       sys,
                          std::shared_ptr<GlobalData>            gd,
                          std::shared_ptr<ExtendedParticleData>  pd,
                          std::vector<std::string> path):sys(sys),gd(gd),pd(pd),path(path){

            integratorsInfo = std::make_shared<InputEntryManager>(sys,path);

            //Load components
            this->loadGroups();
            this->loadSchedule();
            this->loadIntegrators();

            integratorsInfo->checkEntriesUsed();

        }

        IntegratorManager(std::shared_ptr<ExtendedSystem>       sys,
                          std::shared_ptr<GlobalData>            gd,
                          std::shared_ptr<ExtendedParticleData>  pd):IntegratorManager(sys,gd,pd,{"integrator"}){}

        //Add a new interactor to the system
        void addIntegrator(std::shared_ptr<uammd::Integrator> integrator, std::string name, int steps){
            if(integrators.count(name) == 0){
                //Order of added integrators is the number of integrators already added plus one
                stepsInfo info;
                info.name  = name;
                info.order = integratorSteps.size()+1;
                info.steps = steps;

                integratorSteps[name] = info;
                integrators[name]     = integrator;

                System::log<System::MESSAGE>("[IntegratorManager] (%s) Added integrator \"%s\" (manually).",
                                             path.back().c_str(),name.c_str());
                System::log<System::MESSAGE>("[IntegratorManager] (%s) Integrator '%s' will run for %llu steps with order %u.",
                                             path.back().c_str(),name.c_str(),steps,integratorSteps[name].order);

            } else {
                System::log<System::CRITICAL>("[IntegratorManager] (%s) Error adding integrator manually,"
                                              "integrator \"%s\" has already been added.",path.back().c_str(),name.c_str());
            }
        }

        std::map<std::string,std::shared_ptr<uammd::Integrator>>& getIntegrators(){
            return integrators;
        }

        std::shared_ptr<uammd::Integrator> getIntegrator(std::string name){
            if(integrators.count(name) == 0){
                System::log<System::CRITICAL>("[IntegratorManager] (%s) Requested integrator \"%s\" has not been added.",
                                              path.back().c_str(),name.c_str());
            }
            return integrators[name];
        }

        std::vector<stepsInfo> getSortedIntegratorSteps(){
            std::vector<stepsInfo> sortedIntegratorSteps;

            //Load integratror steps into sortedIntegratorSteps
            for(auto& integrator : integratorSteps){
                sortedIntegratorSteps.push_back(integrator.second);
            }

            //Sort sortedIntegratorSteps by order
            std::sort(sortedIntegratorSteps.begin(),sortedIntegratorSteps.end(),
                      [](stepsInfo a, stepsInfo b){return a.order < b.order;});

            //Print information about integrator steps
            for(uint i = 0; i < sortedIntegratorSteps.size(); i++){
                System::log<System::DEBUG>("[IntegratorManager] (%s) Integrator '%s' will run for %llu steps with order %u. (Vector index %u)",
                                           path.back().c_str(),sortedIntegratorSteps[i].name.c_str(),sortedIntegratorSteps[i].steps,
                                           sortedIntegratorSteps[i].order,i);
            }

            //Check if the order of the integrators is correct
            //The first order is 1, which has to be present.
            //The remaining orders must be consecutive.
            int order = 1;
            for(auto& integrator : sortedIntegratorSteps){
                if(integrator.order != order){
                    System::log<System::CRITICAL>("[IntegratorManager] (%s) Integrator '%s' has order %u, but order %u is expected.",
                                                  path.back().c_str(),integrator.name.c_str(),integrator.order,order);
                }
                order++;
            }

            return sortedIntegratorSteps;
        }
};

}}
