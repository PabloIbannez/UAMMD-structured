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
    struct LennardJones_{

        struct ComputationalData{
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

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

            real  r2 = dot(rij, rij);

            real e = LennardJonesType::energy(rij,r2,epsilon,sigma);

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

            real  r2 = dot(rij, rij);

            real3 f = LennardJonesType::force(rij,r2,epsilon,sigma);

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

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;

            const real r2 = dot(rij, rij);

            if        (currentParticleIndex == index_i){
	      H = LennardJonesType::hessian(rij, r2, epsilon, sigma);
            } else if (currentParticleIndex == index_j){
	      H = LennardJonesType::hessian(-rij, r2, epsilon, sigma);
            }

            return H;
        }

    };

    template<class LennardJonesType>
    struct LennardJonesCommon_epsilon_{

        struct ComputationalData: public LennardJones_<LennardJonesType>::ComputationalData{
            real epsilon;
        };

        //Potential parameters

        struct StorageData: public LennardJones_<LennardJonesType>::StorageData{
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
            static_cast<typename LennardJones_<LennardJonesType>::ComputationalData&>(computational) =
            LennardJones_<LennardJonesType>::getComputationalData(gd,pg,storage,computables,st);

            computational.epsilon = storage.epsilon;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.epsilon = data.getParameter<real>("epsilon");

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

            typename LennardJones_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJones_<LennardJonesType>::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            typename LennardJones_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJones_<LennardJonesType>::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            typename LennardJones_<LennardJonesType>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = bondParam.sigma;

            return LennardJones_<LennardJonesType>::hessian(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    struct LennardJonesGaussian_{

        struct ComputationalData{
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real epsilon;
            real sigma;
            real D;
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

            param.epsilon = bondParametersMap.at("epsilon");
            param.sigma   = bondParametersMap.at("sigma");
            param.D       = bondParametersMap.at("D");

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
            const real D       = bondParam.D;

            real r2 = dot(rij, rij);

            real e = BasicPotentials::ModifiedLennardJones::Gaussian::energy(rij,r2,epsilon,sigma,D);

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
            const real D       = bondParam.D;

            real  r2 = dot(rij, rij);

            real3 f = BasicPotentials::ModifiedLennardJones::Gaussian::force(rij,r2,epsilon,sigma,D);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;

        }


    };

    struct LennardJonesGaussianCommon_epsilon_D_ {

        struct ComputationalData: public LennardJonesGaussian_::ComputationalData{
            real epsilon;
            real D;
        };

        //Potential parameters

        struct StorageData: public LennardJonesGaussian_::StorageData{
            real epsilon;
            real D;
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
            static_cast<LennardJonesGaussian_::ComputationalData&>(computational) =
            LennardJonesGaussian_::getComputationalData(gd,pg,storage,computables,st);

            computational.epsilon = storage.epsilon;
            computational.D       = storage.D;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.epsilon = data.getParameter<real>("epsilon");
            storage.D       = data.getParameter<real>("D");

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

            LennardJonesGaussian_::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.D       = computational.D;
            bP.sigma   = bondParam.sigma;

            return LennardJonesGaussian_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            LennardJonesGaussian_::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.D       = computational.D;
            bP.sigma   = bondParam.sigma;

            return LennardJonesGaussian_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    using LennardJonesType1               = Bond2Hessian_<LennardJones_<BasicPotentials::LennardJones::Type1>>;
    using LennardJonesType1Common_epsilon = Bond2Hessian_<LennardJonesCommon_epsilon_<BasicPotentials::LennardJones::Type1>>;

    using LennardJonesType2               = Bond2Hessian_<LennardJones_<BasicPotentials::LennardJones::Type2>>;
    using LennardJonesType2Common_epsilon = Bond2Hessian_<LennardJonesCommon_epsilon_<BasicPotentials::LennardJones::Type2>>;

    using LennardJonesType3               = Bond2Hessian_<LennardJones_<BasicPotentials::LennardJones::Type3>>;
    using LennardJonesType3Common_epsilon = Bond2Hessian_<LennardJonesCommon_epsilon_<BasicPotentials::LennardJones::Type3>>;

    using LennardJonesKaranicolasBrooks               = Bond2_<LennardJones_<BasicPotentials::LennardJones::KaranicolasBrooks>>;
    using LennardJonesKaranicolasBrooksCommon_epsilon = Bond2_<LennardJonesCommon_epsilon_<BasicPotentials::LennardJones::KaranicolasBrooks>>;

    using LennardJonesGaussian                 = Bond2_<LennardJonesGaussian_>;
    using LennardJonesGaussianCommon_epsilon_D = Bond2_<LennardJonesGaussianCommon_epsilon_D_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType1,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType1>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType1Common_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType1Common_epsilon>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType2,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType2>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType2Common_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType2Common_epsilon>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType3,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType3>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesType3Common_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesType3Common_epsilon>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesKaranicolasBrooks,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesKaranicolasBrooks>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesKaranicolasBrooksCommon_epsilon,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesKaranicolasBrooksCommon_epsilon>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesGaussian,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesGaussian>
)

REGISTER_BOND_INTERACTOR(
    Bond2,LennardJonesGaussianCommon_epsilon_D,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::LennardJonesGaussianCommon_epsilon_D>
)
