#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Set/SetInteractor.cuh"
#include "Interactor/Set/Set1/Set1.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Harmonic.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set1{

    struct HarmonicAnisotropic_{

        struct ComputationalData{
            real* mass;
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct SetParameters{
            real3 K;
            real3 r0;
            real3 pos;
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

            param.K   = bondParametersMap.at("K");
            param.r0  = bondParametersMap.at("r0");
            param.pos = bondParametersMap.at("position");

            return param;
        }

        static inline __device__ real energy(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 ref = setParam.pos;
            const real3 rij = computational.box.apply_pbc(ref-centerOfMass);

            const real3 K   = setParam.K;
            const real3 r0  = setParam.r0;

            return BasicPotentials::HarmonicAnisotropic::energy(rij,K,r0);
        }

        static inline __device__ real4 force(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 ref = setParam.pos;
            const real3 rij = computational.box.apply_pbc(ref-centerOfMass);

            const real3 K   = setParam.K;
            const real3 r0  = setParam.r0;

            return make_real4(BasicPotentials::HarmonicAnisotropic::force(rij,K,r0));
        }

    };

    struct HarmonicCommon_K_r0_{


        struct ComputationalData: public HarmonicAnisotropic_::ComputationalData{
            real K;
            real r0;
        };

        //Potential parameters

        struct StorageData: public HarmonicAnisotropic_::StorageData{
            real K;
            real r0;
        };

        struct SetParameters{
            real3 pos;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ExtendedParticleData> pd,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;
            static_cast<HarmonicAnisotropic_::ComputationalData&>(computational) =
            HarmonicAnisotropic_::getComputationalData(gd,pd,storage,comp,st);

            computational.K  = storage.K;
            computational.r0 = storage.r0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ExtendedParticleData> pd,
                                                   DataEntry& data){

            StorageData storage;

            storage.K  = data.getParameter<real>("K");
            storage.r0 = data.getParameter<real>("r0");

            return storage;
        }

        //Set parameters reader

        template<typename T>
        static __host__ SetParameters processSetParameters(std::shared_ptr<GlobalData> gd,
                                                           std::map<std::string,T>& bondParametersMap){

            SetParameters param;

            param.pos = bondParametersMap.at("position");

            return param;
        }

        static inline __device__ real energy(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const typename HarmonicAnisotropic_::SetParameters sP = {{computational.K ,computational.K ,computational.K},
                                                                     {computational.r0,computational.r0,computational.r0},
                                                                      setParam.pos};

            return HarmonicAnisotropic_::energy(centerOfMass,computational,sP);
        }

        static inline __device__ real4 force(const real3& centerOfMass,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const typename HarmonicAnisotropic_::SetParameters sP = {{computational.K ,computational.K ,computational.K},
                                                                     {computational.r0,computational.r0,computational.r0},
                                                                      setParam.pos};

            return HarmonicAnisotropic_::force(centerOfMass,computational,sP);
        }

    };


    using FixedHarmonicAnisotropicCenterOfMass            = CenterOfMass_<HarmonicAnisotropic_>;
    using FixedHarmonicAnisotropicCommon_K_r0CenterOfMass = CenterOfMass_<HarmonicCommon_K_r0_>;

}}}}

REGISTER_SET_INTERACTOR(
    Set1,FixedHarmonicAnisotropicCenterOfMass,
    uammd::structured::Interactor::SetInteractor<uammd::structured::Potentials::Set1::FixedHarmonicAnisotropicCenterOfMass>
)

REGISTER_SET_INTERACTOR(
    Set1,FixedHarmonicAnisotropicCommon_K_r0CenterOfMass,
    uammd::structured::Interactor::SetInteractor<uammd::structured::Potentials::Set1::FixedHarmonicAnisotropicCommon_K_r0CenterOfMass>
)
