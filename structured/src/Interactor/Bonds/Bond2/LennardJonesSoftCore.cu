#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    template<class LennardJonesType>
    struct LennardJonesSoftCore_{

        struct ComputationalData{
            real4* pos;
            Box box;
            real lambda;
            real alpha;
            int  n;
        };

        //Potential parameters

        struct StorageData{
            real alpha;
            int  n;
        };

        struct BondParameters{

            real epsilon;
            real sigma;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box    = gd->getEnsemble()->getBox();
            computational.lambda = gd->getEnsemble()->getLambda();
            computational.alpha  = storage.alpha;
            computational.n      = storage.n;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.alpha = data.getParameter<real>("alpha");
            storage.n     = data.getParameter<int>("n",2);

            System::log<System::MESSAGE>("[LambdaFixedHarmonic] alpha: %f", storage.alpha);
            System::log<System::MESSAGE>("[LambdaFixedHarmonic] n: %i",storage.n);

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.epsilon = bondParametersMap.at("epsilon");
            param.sigma   = bondParametersMap.at("sigma");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            real3 rij = computational.box.apply_pbc(posj-posi);

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;
            const real lambda  = computational.lambda;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            real  r2 = dot(rij, rij);

            real e = LennardJonesType::energy(rij,r2,epsilon,sigma,lambda,alpha,n);

            return e;
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            real3 rij = computational.box.apply_pbc(posj-posi);

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;
            const real lambda  = computational.lambda;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            real  r2 = dot(rij, rij);

            real3 f = LennardJonesType::force(rij,r2,epsilon,sigma,lambda,alpha,n);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;

        }

        static inline __device__ real lambdaDerivative(int index_i, int index_j,
                                                       int currentParticleIndex,
                                                       const ComputationalData &computational,
                                                       const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            real3 rij = computational.box.apply_pbc(posj-posi);

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;
            const real lambda  = computational.lambda;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            real  r2 = dot(rij, rij);

            real ld = LennardJonesType::lambdaDerivative(rij,r2,epsilon,sigma,lambda,alpha,n);

            return ld;
        }
    };

    template<class LennardJonesType>
    struct LennardJonesSoftCoreCommon_epsilon_{

        struct ComputationalData: public LennardJonesSoftCore_<LennardJonesType>::ComputationalData{
            real epsilon;
        };

        //Potential parameters

        struct StorageData: public LennardJonesSoftCore_<LennardJonesType>::StorageData{
            real epsilon;
        };

        struct BondParameters{
            real sigma;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;
            static_cast<typename LennardJonesSoftCore_<LennardJonesType>::ComputationalData&>(computational) =
            LennardJonesSoftCore_<LennardJonesType>::getComputationalData(gd,pg,storage,computables,st);

            computational.epsilon = storage.epsilon;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;
            static_cast<typename LennardJonesSoftCore_<LennardJonesType>::StorageData&>(storage) =
            LennardJonesSoftCore_<LennardJonesType>::getStorageData(gd,pg,data);

            storage.epsilon = data.getParameter<real>("epsilon");

            System::log<System::MESSAGE>("[LambdaFixedHarmonicCommon_epsilon] epsilon: %f", storage.epsilon);

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.sigma   = bondParametersMap.at("sigma");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            typename LennardJonesSoftCore_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJonesSoftCore_<LennardJonesType>::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            typename LennardJonesSoftCore_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJonesSoftCore_<LennardJonesType>::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

        static inline __device__ real lambdaDerivative(int index_i, int index_j,
                                                       int currentParticleIndex,
                                                       const ComputationalData &computational,
                                                       const BondParameters &bondParam){

            typename LennardJonesSoftCore_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJonesSoftCore_<LennardJonesType>::lambdaDerivative(index_i,index_j,currentParticleIndex,computational,bP);
        }

    };

    using LennardJonesSoftCoreType1 = Bond2Lambda_<LennardJonesSoftCore_<BasicPotentials::LennardJones::SoftCoreType1>>;
    using LennardJonesSoftCoreType2 = Bond2Lambda_<LennardJonesSoftCore_<BasicPotentials::LennardJones::SoftCoreType2>>;

    using LennardJonesSoftCoreType1Common_epsilon = Bond2Lambda_<LennardJonesSoftCoreCommon_epsilon_<BasicPotentials::LennardJones::SoftCoreType1>>;
    using LennardJonesSoftCoreType2Common_epsilon = Bond2Lambda_<LennardJonesSoftCoreCommon_epsilon_<BasicPotentials::LennardJones::SoftCoreType2>>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesSoftCoreType1,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesSoftCoreType1>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesSoftCoreType2,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesSoftCoreType2>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesSoftCoreType1Common_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesSoftCoreType1Common_epsilon>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesSoftCoreType2Common_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesSoftCoreType2Common_epsilon>
)
