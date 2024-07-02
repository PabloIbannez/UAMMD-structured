#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "External.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

#include "ExternalTabulated.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct PolarizationTabulated_: public ExternalTabulated_{

        using StorageData = ExternalTabulated_::StorageData;

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data) {
            return ExternalTabulated_::getStorageData(gd, pg, data);
        };

        struct ComputationalData : public ExternalTabulated_::ComputationalData {
            real* polarizability;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            static_cast<ExternalTabulated_::ComputationalData&>(computational) =
            ExternalTabulated_::getComputationalData(gd, pg, storage, computables, st);

            // Set polarizability

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.polarizability = pd->getPolarizability(access::location::gpu, access::mode::read).raw();

            return computational;
        };

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
                real e = ExternalTabulated_::energy(index_i, computational);
                return e*computational.polarizability[index_i];
        }

        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
                real3 f = ExternalTabulated_::force(index_i, computational);
                return f*computational.polarizability[index_i];
        }

    };

    using PolarizationTabulated   = External_<PolarizationTabulated_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    External,PolarizationTabulated,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::PolarizationTabulated>
)
