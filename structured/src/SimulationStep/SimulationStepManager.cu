#include "SimulationStep/SimulationStepManager.cuh"

namespace uammd{
namespace structured{

    void SimulationStepManager::loadGroups(){
        groups = GroupUtils::loadGroupsListFromInputEntries(sys,gd,pd,simulationStepsInfo);
    }

    void SimulationStepManager::loadSimulationSteps(){

        for(auto& entry : simulationStepsInfo->getEntriesInfo()){
            if(SimulationStepLoader::isSimulationStepAvailable(sys,entry.second.path)){
                std::shared_ptr<SimulationStep::SimulationStepBase> simStep = SimulationStepLoader::loadSimulationStep(sys,groups,integrator,ff,entry.second.path);

                if(simulationSteps.count(entry.second.name) == 0){

                    simulationSteps[entry.second.name] = simStep;
                    entry.second.used = true;
                }
                else{
                    System::log<System::CRITICAL>("[SimulationStepManager] (%s) Error loading simulation step,"
                                                  "simulation step \"%s\" has already been loaded.",path.back().c_str(),entry.second.name.c_str());
                }

            }
        }

        //Print information about loaded simulation steps
        for(auto& simS : simulationSteps){
            System::log<System::MESSAGE>("[SimulationStepManager] (%s) Simulation step '%s' loaded.",
                                         path.back().c_str(),simS.first.c_str());
        }
    }

    SimulationStepManager::SimulationStepManager(std::shared_ptr<IntegratorManager>  integrator,
                                                 std::shared_ptr<ForceField>                 ff,
                                                 std::vector<std::string> path):integrator(integrator),ff(ff),path(path){

        topology = ff->getTopology();

        sys = topology->getSystem();
        gd  = topology->getGlobalData();
        pd  = topology->getParticleData();


        simulationStepsInfo = std::make_shared<InputEntryManager>(sys,path);

        //Load components
        this->loadGroups();
        this->loadSimulationSteps();

        simulationStepsInfo->checkEntriesUsed();

    }

    SimulationStepManager::SimulationStepManager(std::shared_ptr<IntegratorManager>  integrator,
                                                 std::shared_ptr<ForceField>                 ff):
    SimulationStepManager(integrator,ff,{"simulationStep"}){}

    //Add a new interactor to the system
    void SimulationStepManager::addSimulationStep(std::shared_ptr<SimulationStep::SimulationStepBase> simStep,
                                                  std::string name){
        if(simulationSteps.count(name) == 0){

            simulationSteps[name]     = simStep;

            System::log<System::MESSAGE>("[SimulationStepManager] (%s) Added simulation step \"%s\" (manually).",
                                         path.back().c_str(),name.c_str());
        } else {
            System::log<System::CRITICAL>("[SimulationStepManager] (%s) Error adding simulation step manually,"
                                          "simulation step \"%s\" has already been added.",path.back().c_str(),name.c_str());
        }
    }

}}
