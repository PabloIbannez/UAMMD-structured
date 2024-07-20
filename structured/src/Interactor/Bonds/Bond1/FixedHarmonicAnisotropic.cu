#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond1/Bond1.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Harmonic.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond1{

    struct FixedHarmonicAnisotropic_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{
            real3 K;
            real3 r0;
            real3 pos;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox();

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

            param.K    = bondParametersMap.at("K");
            param.r0   = bondParametersMap.at("r0");
            param.pos  = bondParametersMap.at("position");

            return param;
        }

        //Energy and force definition

        static inline __device__ real3 force(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

            return BasicPotentials::HarmonicAnisotropic::force(rij,K,r0);
        }

        static inline __device__ real energy(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

            const real e = BasicPotentials::HarmonicAnisotropic::energy(rij,K,r0);

            return e;
        }

      static inline __device__ tensor3 hessian(int index_i,
					       int currentParticleIndex,
					       const ComputationalData &computational,
					       const BondParameters &bondParam){

	tensor3 H = tensor3();
	H.xx = bondParam.K.x; H.yy = bondParam.K.y; H.zz = bondParam.K.z;
	return H;
        }

    };

    using FixedHarmonicAnisotropic = Bond1Hessian_<FixedHarmonicAnisotropic_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond1,FixedHarmonicAnisotropic,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond1::FixedHarmonicAnisotropic>
)
