#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationUtils{

class SortStep: public SimulationStepBase {

    public:

        SortStep(std::shared_ptr<ParticleGroup>             pg,
                 std::shared_ptr<IntegratorManager> integrator,
                 std::shared_ptr<ForceField>                ff,
                 DataEntry& data,
                 std::string name):SimulationStepBase(pg,integrator,ff,data,name){}

        void init(cudaStream_t st) override {}

        void applyStep(ullint step, cudaStream_t st) override {
            this->pd->sortParticles(st);
        }

};

}}}}

REGISTER_SIMULATION_STEP(
    UtilsStep,SortStep,
    uammd::structured::SimulationStep::SimulationUtils::SortStep
)
