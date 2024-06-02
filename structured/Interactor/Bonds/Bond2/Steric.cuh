#ifndef __STERIC_BOND2__
#define __STERIC_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    template<int power>
    struct Steric_{

        struct ComputationalData{
            real4* pos;
            Box    box;
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
                                                               const Computables& computables){

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

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;

            const real r2 = dot(rij, rij);

            return BasicPotentials::Steric::Steric::energy<power>(rij,r2,epsilon,sigma);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real epsilon = bondParam.epsilon;
            const real sigma   = bondParam.sigma;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::Steric::Steric::force<power>(rij,r2,epsilon,sigma);

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
	  H = BasicPotentials::Steric::Steric::hessian<power>(rij, r2, epsilon, sigma);
	} else if (currentParticleIndex == index_j){
	  H = BasicPotentials::Steric::Steric::hessian<power>(-rij, r2, epsilon, sigma);
	}

	return H;
        }
    };

    template<int power>
    struct StericCommon_epsilon_sigma_{

        struct ComputationalData: public Steric_<power>::ComputationalData{
            real epsilon;
            real sigma;
        };

        //Potential parameters

        struct StorageData: public Steric_<power>::StorageData{
            real epsilon;
            real sigma;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables){
            ComputationalData computational;
            static_cast<typename Steric_<power>::ComputationalData&>(computational) = Steric_<power>::getComputationalData(gd,pg,storage,computables);

            computational.epsilon = storage.epsilon;
            computational.sigma   = storage.sigma;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.epsilon = data.getParameter<real>("epsilon");
            storage.sigma   = data.getParameter<real>("sigma");

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

            typename Steric_<power>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = computational.sigma;

            return Steric_<power>::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            typename Steric_<power>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = computational.sigma;

            return Steric_<power>::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
					       int currentParticleIndex,
					       const ComputationalData &computational,
					       const BondParameters &bondParam){

            typename Steric_<power>::BondParameters bP;
            bP.epsilon = computational.epsilon;
            bP.sigma   = computational.sigma;

            return Steric_<power>::hessian(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    using Steric6                      = Bond2Hessian_<Steric_<6>>;
    using Steric6Common_epsilon_sigma  = Bond2Hessian_<StericCommon_epsilon_sigma_<6>>;

    using Steric12                     = Bond2Hessian_<Steric_<12>>;
    using Steric12Common_epsilon_sigma = Bond2Hessian_<StericCommon_epsilon_sigma_<12>>;

}}}}

#endif
