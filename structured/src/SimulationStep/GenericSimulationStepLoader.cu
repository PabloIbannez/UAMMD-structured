#include "SimulationStep/GenericSimulationStepLoader.cuh"

namespace uammd{
namespace structured{
namespace SimulationStepLoader{

    bool isSimulationStepAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string simulationStepType    = data.getType();
        std::string simulationStepSubType = data.getSubType();

        return SimulationStep::SimulationStepFactory::getInstance().isSimulationStepRegistered(simulationStepType,simulationStepSubType);

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

        return SimulationStep::SimulationStepFactory::getInstance().createSimulationStep(simulationStepType,simulationStepSubType,
                                                                                         pg,integrator,ff,
                                                                                         data,
                                                                                         path.back());

    }

}}}
