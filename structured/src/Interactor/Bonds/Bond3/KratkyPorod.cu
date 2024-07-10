#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond3/Bond3.cuh"
#include "Interactor/InteractorFactory.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond3{

    struct KratkyPorod_{

        //Potential parameters

        struct ComputationalData{
            Box    box;
            real4* pos;
        };

        struct StorageData{};

        struct BondParameters{
            real K;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.box = gd->getEnsemble()->getBox();
            computational.pos = pd->getPos(access::location::gpu,
                                           access::mode::read).raw();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){
            StorageData storage;
            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.K      = bondParametersMap.at("K");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(const real& ang,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            return bondParam.K*(real(1.0)+cos(ang));
        }

        static inline __device__ real energyDerivate(const real& ang,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){
            return -bondParam.K*sin(ang);
        }

      static inline __device__ real energySecondDerivate(const real& ang,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

	return -bondParam.K*cos(ang);
      }
    };

    struct KratkyPorodCommon_K_{

        //Potential parameters
        struct ComputationalData : public KratkyPorod_::ComputationalData{
            real K;
        };

        struct StorageData : public KratkyPorod_::StorageData{
            real K;
        };

        struct BondParameters{};

        //Storage data reader

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){
            ComputationalData computational;
            static_cast<KratkyPorod_::ComputationalData&>(computational) =
            KratkyPorod_::getComputationalData(gd, pg, storage, computables, st);

            computational.K = storage.K;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K   = data.getParameter<real>("K");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            return param;
        }

        //Angular potential functions

        static inline __device__ real energy(const real& ang,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            KratkyPorod_::BondParameters bP;
            bP.K      = computational.K;

            return KratkyPorod_::energy(ang,computational,bP);
        }

        static inline __device__ real energyDerivate(const real& ang,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            KratkyPorod_::BondParameters bP;
            bP.K      = computational.K;

            return KratkyPorod_::energyDerivate(ang,computational,bP);
	}

      static inline __device__ real energySecondDerivate(const real& ang,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

	KratkyPorod_::BondParameters bP;
	bP.K      = computational.K;

	return KratkyPorod_::energySecondDerivate(ang,computational,bP);
      }
    };

    using KratkyPorod         = Bond3Hessian_<KratkyPorod_>;
    using KratkyPorodCommon_K = Bond3Hessian_<KratkyPorodCommon_K_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond3,KratkyPorod,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond3::KratkyPorod>
)

REGISTER_BOND_INTERACTOR(
    Bond3,KratkyPorodCommon_K,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond3::KratkyPorodCommon_K>
)

