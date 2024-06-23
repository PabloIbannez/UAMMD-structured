#ifndef __SIMULATION_STEP_LOADER__
#define __SIMULATION_STEP_LOADER__
namespace uammd{
namespace structured{
namespace SimulationStepLoader{

    bool isSimulationStepAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();
        
        return false;

    }

    
    std::shared_ptr<SimulationStep::SimulationStepBase>
    loadSimulationStep(std::shared_ptr<ExtendedSystem> sys,
                       std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField> ff,
                       std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();

        std::shared_ptr<SimulationStep::SimulationStepBase> simulationStep;
        bool found = false;
        

        if(not found){
            System::log<System::CRITICAL>("[SimulationStepLoader] (%s) Could not find simulationStep %s::%s",
                                            path.back().c_str(),simulationStepType.c_str(),simulationStepSubType.c_str());
        }

        return simulationStep;

    }

    }}}
#endif
