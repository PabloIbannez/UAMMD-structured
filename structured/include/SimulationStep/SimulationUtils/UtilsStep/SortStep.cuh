#ifndef __SORT_STEP__
#define __SORT_STEP__

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

#endif
