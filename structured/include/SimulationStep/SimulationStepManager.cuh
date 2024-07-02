#pragma once

#include "Integrator/IntegratorManager.cuh"
#include "ForceFields/ForceFields.cuh"

#include "Topology/Topology.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/GenericSimulationStepLoader.cuh"

namespace uammd{
namespace structured{

class SimulationStepManager{

    private:

        std::shared_ptr<IntegratorManager>  integrator;
        std::shared_ptr<ForceField>                 ff;

        std::shared_ptr<Topology>        topology;

        std::shared_ptr<ExtendedSystem>       sys;
        std::shared_ptr<GlobalData>            gd;
        std::shared_ptr<ExtendedParticleData>  pd;

        std::vector<std::string> path;

        //////////////////////////////////////////

        std::shared_ptr<InputEntryManager>   simulationStepsInfo;

        //////////////////////////////////////////

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;
        std::map<std::string,std::shared_ptr<SimulationStep::SimulationStepBase>> simulationSteps;

        //////////////////////////////////////////

        void loadGroups();
        void loadSimulationSteps();

    public:

        SimulationStepManager(std::shared_ptr<IntegratorManager>  integrator,
                              std::shared_ptr<ForceField>                 ff,
                              std::vector<std::string> path);

        SimulationStepManager(std::shared_ptr<IntegratorManager>  integrator,
                              std::shared_ptr<ForceField>                 ff);

        //Add a new interactor to the system
        void addSimulationStep(std::shared_ptr<SimulationStep::SimulationStepBase> simStep, std::string name);

        std::map<std::string,std::shared_ptr<SimulationStep::SimulationStepBase>>& getSimulationSteps(){
            return simulationSteps;
        }

        std::shared_ptr<SimulationStep::SimulationStepBase> getSimulationStep(std::string name){
            if(simulationSteps.count(name) == 0){
                System::log<System::CRITICAL>("[SimulationStepManager] (%s) Requested simulation step \"%s\" has not been added.",
                                              path.back().c_str(),name.c_str());
            }
            return simulationSteps[name];
        }
};

}}
