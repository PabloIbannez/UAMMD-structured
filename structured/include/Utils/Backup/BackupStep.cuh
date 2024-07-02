#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"
#include "Integrator/IntegratorManager.cuh"

#include "ForceFields/ForceFields.cuh"

#include "SimulationStep/SimulationStep.cuh"

#include <memory>
#include <string>

namespace uammd{
namespace structured{
namespace Backup{

class BackupStep: public SimulationStep::SimulationStepBase {

    private:

        std::shared_ptr<ExtendedSystem::InputType> input;

        bool isSimulationStateCorrect();

        void writeBackup(std::string backupFileName);

    public:

        struct Parameters : public SimulationStep::SimulationStepBase::Parameters {};

        BackupStep(std::shared_ptr<ParticleGroup>              pg,
                   std::shared_ptr<IntegratorManager>  integrator,
                   std::shared_ptr<ForceField>                 ff,
                   Parameters par,
                   std::string name);

        ~BackupStep();

        void init(cudaStream_t st) override;

        void applyStep(ullint step, cudaStream_t st) override;

};

}}}

