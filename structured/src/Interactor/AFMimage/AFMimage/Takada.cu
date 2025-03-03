#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/AFMimage/AFMimageInteractor.cuh"
#include "Interactor/AFMimage/AFMimage/AFMimage.cuh"
#include "Interactor/InteractorFactory.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace AFMimage{

    struct Takada_{

        struct ComputationalData {};

        struct StorageData {};

        //Computational data getter
        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                      std::shared_ptr<ParticleGroup> pg,
                                                      const StorageData&  storage,
                                                      const Computables& comp,
                                                      const cudaStream_t& st){
            ComputationalData computational;
            return computational;
        }

        //Storage data reader
        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          DataEntry& data){

            StorageData storage;
            return storage;
        }

        //Energy and force definition

        static inline __device__ real energy(int index,
                                             int tipIndex,
                                             const ComputationalData &computational){
            const real e = real(0.0);
            return e;
        }

        static inline __device__ real3 force(int index,
                                             int tipIndex,
                                             const ComputationalData &computational){
            real3 f = make_real3(0.0);
            return f;

        }

    };

    using Takada = AFMimage_<Takada_>;

}}}}

REGISTER_AFM_IMAGE_INTERACTOR(
    AFMimage,Takada,
    uammd::structured::Interactor::AFMimageInteractor<uammd::structured::Potentials::AFMimage::Takada>
)
