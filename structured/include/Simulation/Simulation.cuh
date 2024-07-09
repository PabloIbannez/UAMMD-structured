#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"

#include "Topology/Topology.cuh"
#include "Integrator/IntegratorManager.cuh"

#include "ForceFields/ForceFields.cuh"

#include "SimulationStep/SimulationStepManager.cuh"

#include "Utils/Backup/BackupStep.cuh"

#include<unistd.h>
#include<sys/wait.h>

namespace uammd{
namespace structured{

class Simulation{

    private:

        std::shared_ptr<ExtendedSystem>      sys;

        std::shared_ptr<GlobalData>           gd;
        std::shared_ptr<ExtendedParticleData> pd;

        /////////////////////////////////////////

        std::shared_ptr<Topology>           topology;
        std::shared_ptr<IntegratorManager>  integrators;

        std::shared_ptr<ForceField>         ff;

        std::shared_ptr<SimulationStepManager> simulationSteps;

    public:

        Simulation(std::shared_ptr<ExtendedSystem> sys);

        //Getters

        std::shared_ptr<ExtendedSystem>       getSystem(){ return this->sys;}
        std::shared_ptr<GlobalData>           getGlobalData(){ return this->gd;}
        std::shared_ptr<ExtendedParticleData> getParticleData(){return this->pd;}

        std::shared_ptr<Topology>   getTopology(){return this->topology;}
        std::shared_ptr<ForceField> getForceField(){return this->ff;}

        std::shared_ptr<IntegratorManager>   getIntegrators(){return this->integrators;}

        std::shared_ptr<SimulationStepManager> getSimulationSteps(){return this->simulationSteps;}

        //Run simulation

        int run();
};

template<class inTyp>
void startSelfStartingSimulation(const inTyp& in);

void startSelfStartingSimulationFromFile(std::string inputFilePath);

void startSelfStartingSimulationFromInput(const typename Input::Input::DataType& in);

}}
