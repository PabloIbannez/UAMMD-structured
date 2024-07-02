#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Bond3.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond3{

    struct BestChenHummerAngular_{

        private:

            static constexpr real gamma = real(0.1);

            static constexpr real k_alpha = real(106.4);
            static constexpr real k_beta  = real(26.3);

            static constexpr real theta_alpha = real(1.60);
            static constexpr real theta_beta  = real(2.27);

            static constexpr real epsilon_alpha = real(4.3);

        public:

            //Potential parameters

            struct ComputationalData{
                Box    box;
                real4* pos;
            };

            struct StorageData{};

            struct BondParameters{};

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

            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                       std::shared_ptr<ParticleGroup>        pg,
                                                       DataEntry& data){

                StorageData storage;

                //Check if the used units are KcalMol_A
                std::string unitsName = gd->getUnits()->getSubType();
                if(unitsName != "KcalMol_A"){
                    System::log<System::CRITICAL>("[BestChenHummerAngular] The potential is only defined for KcalMol_A units."
                                                  " But the used units are %s", unitsName.c_str());
                }

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

                real adiff_alpha2 = ang-theta_alpha;
                     adiff_alpha2 = adiff_alpha2*adiff_alpha2;
                const real exp_alpha = exp(-gamma*(k_alpha*adiff_alpha2+epsilon_alpha));

                real adiff_beta2 = ang-theta_beta;
                     adiff_beta2 = adiff_beta2*adiff_beta2;
                const real exp_beta = exp(-gamma*k_beta*adiff_beta2);

                return -real(1.0/gamma)*log(exp_alpha+exp_beta);
            }


            static inline __device__ real energyDerivate(const real& ang,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

                const real adiff_alpha  = ang-theta_alpha;
                const real adiff_alpha2 = adiff_alpha*adiff_alpha;
                const real exp_alpha = exp(-gamma*(k_alpha*adiff_alpha2+epsilon_alpha));

                const real adiff_beta  = ang-theta_beta;
                const real adiff_beta2 = adiff_beta*adiff_beta;
                const real exp_beta = exp(-gamma*k_beta*adiff_beta2);

                return real(2.0)*(k_alpha*adiff_alpha*exp_alpha+k_beta*adiff_beta*exp_beta)/(exp_alpha+exp_beta);
            }

      static inline __device__ real energySecondDerivate(const real& ang,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

                const real adiff_alpha  = ang-theta_alpha;
                const real adiff_alpha2 = adiff_alpha*adiff_alpha;
                const real exp_alpha = exp(-gamma*(k_alpha*adiff_alpha2+epsilon_alpha));

                const real adiff_beta  = ang-theta_beta;
                const real adiff_beta2 = adiff_beta*adiff_beta;
                const real exp_beta = exp(-gamma*k_beta*adiff_beta2);

		const real deriv_exp_alpha = -real(2.0)*gamma*k_alpha*adiff_alpha*exp_alpha;
		const real deriv_exp_beta  = -real(2.0)*gamma*k_beta*adiff_beta*exp_beta;

		const real term_alpha = real(2.0)*gamma*k_alpha*(deriv_exp_alpha*adiff_alpha + exp_alpha);
		const real term_beta  = real(2.0)*gamma*k_beta*(deriv_exp_beta*adiff_beta + exp_beta);

		const real term_alpha_beta = (deriv_exp_alpha + deriv_exp_beta)*(deriv_exp_alpha + deriv_exp_beta)/(exp_alpha+exp_beta);
		return (term_alpha + term_beta + term_alpha_beta)/(gamma * (exp_alpha + exp_beta));
    }

    };

    using BestChenHummerAngular = Bond3Hessian_<BestChenHummerAngular_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond3,BestChenHummerAngular,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond3::BestChenHummerAngular>
)
