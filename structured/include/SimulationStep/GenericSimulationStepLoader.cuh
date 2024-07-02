#pragma once

#include"SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStepLoader{

    bool isSimulationStepAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path);


    std::shared_ptr<SimulationStep::SimulationStepBase>
    loadSimulationStep(std::shared_ptr<ExtendedSystem> sys,
                       std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField> ff,
                       std::vector<std::string>       path);

}}}
