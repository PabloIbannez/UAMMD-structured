#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Harmonic.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct Harmonic_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real K;
            real r0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
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

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
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

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            const real r2 = dot(rij, rij);

            const real e = BasicPotentials::Harmonic::energy(rij,r2,K,r0);

            return e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::Harmonic::force(rij,r2,K,r0);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }

        static inline __device__ tensor3 hessian(int index_i, int index_j,
                                                 int currentParticleIndex,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam){

            tensor3 H = tensor3(0.0);

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);
	    const real r2   = dot(rij, rij);

            if        (currentParticleIndex == index_i){
	      H = BasicPotentials::Harmonic::hessian(rij,r2,bondParam.K,bondParam.r0);
            } else if (currentParticleIndex == index_j){
	      H = BasicPotentials::Harmonic::hessian(-rij,r2,bondParam.K,bondParam.r0);
            }

            return H;
        }

    };

    struct HarmonicCommon_K_{

        struct ComputationalData : public Harmonic_::ComputationalData{
            real K;
        };

        //Potential parameters

        struct StorageData : public Harmonic_::StorageData{
            real K;
        };

        struct BondParameters{
            real r0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;
            static_cast<Harmonic_::ComputationalData&>(computational)
            = Harmonic_::getComputationalData(gd,pg,storage,computables,st);

            computational.K = storage.K;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K = data.getParameter<real>("K");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.r0  = bondParametersMap.at("r0");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = bondParam.r0;

            return Harmonic_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = bondParam.r0;

            return Harmonic_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

        static inline __device__ tensor3 hessian(int index_i, int index_j,
                                                 int currentParticleIndex,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = bondParam.r0;

            return Harmonic_::hessian(index_i,index_j,currentParticleIndex,computational,bP);
        }

    };

    struct HarmonicCommon_K_r0_{

        struct ComputationalData : public Harmonic_::ComputationalData{
            real K;
            real r0;
        };

        //Potential parameters

        struct StorageData : public Harmonic_::StorageData{
            real K;
            real r0;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){
            ComputationalData computational;
            static_cast<Harmonic_::ComputationalData&>(computational)
            = Harmonic_::getComputationalData(gd,pg,storage,computables,st);

            computational.K  = storage.K;
            computational.r0 = storage.r0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.r0  = data.getParameter<real>("r0");
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

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = computational.r0;

            return Harmonic_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = computational.r0;

            return Harmonic_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

        static inline __device__ tensor3 hessian(int index_i, int index_j,
                                                 int currentParticleIndex,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam){

            Harmonic_::BondParameters bP;
            bP.K  = computational.K;
            bP.r0 = computational.r0;

            return Harmonic_::hessian(index_i,index_j,currentParticleIndex,computational,bP);
        }

    };

    using Harmonic            = Bond2Hessian_<Harmonic_>;
    using HarmonicCommon_K    = Bond2Hessian_<HarmonicCommon_K_>;
    using HarmonicCommon_K_r0 = Bond2Hessian_<HarmonicCommon_K_r0_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,Harmonic,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::Harmonic>
)

REGISTER_BOND_INTERACTOR(
    Bond2,HarmonicCommon_K,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::HarmonicCommon_K>
)

REGISTER_BOND_INTERACTOR(
    Bond2,HarmonicCommon_K_r0,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::HarmonicCommon_K_r0>
)
