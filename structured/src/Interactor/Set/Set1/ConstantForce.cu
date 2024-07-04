#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Set/SetInteractor.cuh"
#include "Interactor/Set/Set1/Set1.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set1{

    struct ConstantForce_{

        struct ComputationalData{
            real* mass;
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct SetParameters{
            real3 F;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ExtendedParticleData> pd,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            computational.pos  = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.mass = pd->getMass(access::location::gpu, access::mode::read).raw();

            computational.box  = gd->getEnsemble()->getBox();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ExtendedParticleData> pd,
                                                   DataEntry& data){

            StorageData storage;
            return storage;
        }

        //Set parameters reader

        template<typename T>
        static __host__ SetParameters processSetParameters(std::shared_ptr<GlobalData> gd,
                                                           std::map<std::string,T>& bondParametersMap){

            SetParameters param;

            param.F    = bondParametersMap.at("force");

            return param;
        }

        static inline __device__ real energy(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            return real(0.0);
        }

        static inline __device__ real4 force(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            return make_real4(setParam.F,0.0);
        }

    };

    using ConstantForceOverCenterOfMass = CenterOfMass_<ConstantForce_>;

}}}}

REGISTER_SET_INTERACTOR(
    Set1,ConstantForceOverCenterOfMass,
    uammd::structured::Interactor::SetInteractor<uammd::structured::Potentials::Set1::ConstantForceOverCenterOfMass>
)
