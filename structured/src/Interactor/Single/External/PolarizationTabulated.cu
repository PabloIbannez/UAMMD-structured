#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/External/External.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "ExternalTabulated.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct PolarizationTabulated_: public ExternalTabulated_{

        struct StorageData : public ExternalTabulated_::StorageData {
            real eps;

        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data) {

            StorageData storage;

            static_cast<ExternalTabulated_::StorageData&>(storage) =
                ExternalTabulated_::getStorageData(gd,pg,data);

            // Set eps

            real f = gd->getUnits()->getElectricConversionFactor();
            real dielectricConstant = data.getParameter<real>("dielectricConstant");
            storage.eps = dielectricConstant/(4*M_PI*f);

            return storage;
        };

        struct ComputationalData : public ExternalTabulated_::ComputationalData {
            real* polarizability;
            real eps;

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

            // Set eps

            computational.eps = storage.eps;

            return computational;
        };

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
                real e = ExternalTabulated_::energy(index_i, computational);
                return e*computational.eps*computational.polarizability[index_i];
        }

        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
                real3 f = ExternalTabulated_::force(index_i, computational);
                return f*computational.eps*computational.polarizability[index_i];
        }

    };

    using PolarizationTabulated   = External_<PolarizationTabulated_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    External,PolarizationTabulated,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::PolarizationTabulated>
)
