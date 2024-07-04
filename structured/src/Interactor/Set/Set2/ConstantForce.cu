#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Set/SetInteractor.cuh"
#include "Interactor/Set/Set2/Set2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set2{

    struct ConstantForce_{

        struct ComputationalData{
            real4* pos;
            real* mass;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct SetParameters{
            real F;
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

            computational.box = gd->getEnsemble()->getBox();

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

        static inline __device__ real energy(const real3& centerOfMass1,
                                             const real3& centerOfMass2,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 r12 = computational.box.apply_pbc(centerOfMass2-centerOfMass1);
            const real  r2  = dot(r12,r12);

            return BasicPotentials::ConstantForce::energy(r12,r2,setParam.F);
        }

        static inline __device__ real4 force(const real3& centerOfMass1,
                                             const real3& centerOfMass2,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 r12 = computational.box.apply_pbc(centerOfMass2-centerOfMass1);
            const real  r2  = dot(r12,r12);

            return make_real4(BasicPotentials::ConstantForce::force(r12,r2,setParam.F),real(0.0));
        }

    };

    using ConstantForceBetweenCentersOfMass = CenterOfMass_<ConstantForce_>;

}}}}

REGISTER_SET_INTERACTOR(
    Set2,ConstantForceBetweenCentersOfMass,
    uammd::structured::Interactor::SetInteractor<uammd::structured::Potentials::Set2::ConstantForceBetweenCentersOfMass>
)
