#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"
#include "Interactor/BasicParameters.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    struct SplitLennardJones_{

        using LennardJonesType      = typename BasicPotentials::LennardJones::Type1;

        using ParametersType        = typename BasicParameters::Pairs::LennardJones;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        struct ComputationalData{

            real4* pos;
            Box box;

            ParametersPairsIterator paramPairIterator;

            real cutOffFactor;

            real epsilon_r;
            real epsilon_a;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> ljParam;

            real cutOffFactor;
            real cutOff;

            real epsilon_r;
            real epsilon_a;
        };

        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox();

            computational.paramPairIterator = storage.ljParam->getPairIterator();

            computational.cutOffFactor = storage.cutOffFactor;

            computational.epsilon_r = storage.epsilon_r;
            computational.epsilon_a = storage.epsilon_a;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;
            storage.cutOffFactor  = data.getParameter<real>("cutOffFactor");

            storage.epsilon_r = data.getParameter<real>("epsilon_r");
            storage.epsilon_a = data.getParameter<real>("epsilon_a");

            ///////////////////////////////////////////////////////////

            storage.ljParam = std::make_shared<ParameterPairsHandler>(gd,pg,
                                                                    data);

            ///////////////////////////////////////////////////////////

            auto pairsParam = storage.ljParam->getPairParameters();

            real maxSigma = 0.0;
            for(auto p : pairsParam){
                maxSigma=std::max(maxSigma,p.second.sigma);
            }

            storage.cutOff = maxSigma*storage.cutOffFactor;

            System::log<System::MESSAGE>("[LennardJones] cutOff: %f" ,storage.cutOff);

            return storage;
        }

        static inline __device__ real energy(int index_i, int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon_r = computational.epsilon_r;
            const real epsilon_a = computational.epsilon_a;

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

            const real r2 = dot(rij, rij);

            real e = real(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){

                //r2 > (sigma*2^(1/6))^2
                real near_cut2 = sigma*sigma*real(1.259921);

                real eps;
                real eps_a = epsilon*epsilon_a;
                if(r2 < near_cut2){
                    eps = epsilon_r;
                } else {
                    eps = eps_a;
                }

                e = LennardJonesType::energy(rij,r2,eps,sigma) + eps - eps_a;
            }

            return e;

        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon_r = computational.epsilon_r;
            const real epsilon_a = computational.epsilon_a;

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

            const real r2 = dot(rij, rij);

            real3 f = make_real3(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){

                //r2 > (sigma*2^(1/6))^2
                real near_cut2 = sigma*sigma*real(1.259921);

                real eps;
                real eps_a = epsilon*epsilon_a;
                if(r2 < near_cut2){
                    eps = epsilon_r;
                } else {
                    eps = eps_a;
                }

                f = LennardJonesType::force(rij,r2,eps,sigma);

            }

            return f;
        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
					       const ComputationalData& computational){

	const real4 posi = computational.pos[index_i];
	const real4 posj = computational.pos[index_j];

	const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

	const real epsilon_r = computational.epsilon_r;
	const real epsilon_a = computational.epsilon_a;

	const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
	const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

	const real r2 = dot(rij, rij);

	tensor3 H = tensor3(0.0);

	real cutOff2=sigma*computational.cutOffFactor;
	cutOff2=cutOff2*cutOff2;
	if(r2<=cutOff2){
	  real near_cut2 = sigma*sigma*real(1.259921);

	  real eps;
	  real eps_a = epsilon*epsilon_a;
	  if(r2 < near_cut2){
	    eps = epsilon_r;
	  } else {
	    eps = eps_a;
	  }
	  H = LennardJonesType::hessian(rij,r2,eps,sigma);
	}
	return H;
      }


    };

    using SplitLennardJones = NonBondedHessian_<SplitLennardJones_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,SplitLennardJones,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::SplitLennardJones>
)
