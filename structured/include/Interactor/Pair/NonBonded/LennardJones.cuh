#ifndef __LENNARD_JONES_POTENTIAL__
#define __LENNARD_JONES_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template<class LennardJonesType>
    struct LennardJones_{

        using ParametersType        = typename BasicParameters::Pairs::LennardJones;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        struct ComputationalData{

            real4* pos;
            Box box;

            ParametersPairsIterator paramPairIterator;

            real cutOffFactor;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> ljParam;

            real cutOffFactor;
            real cutOff;
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

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;
            storage.cutOffFactor  = data.getParameter<real>("cutOffFactor");

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

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

            const real r2 = dot(rij, rij);

            real e = real(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                e = LennardJonesType::energy(rij,r2,epsilon,sigma);
            }

            return e;

        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

            const real r2 = dot(rij, rij);

            real3 f = make_real3(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                f = LennardJonesType::force(rij,r2,epsilon,sigma);
            }

            return f;
        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
					       const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;

            const real r2 = dot(rij, rij);

            tensor3 H = tensor3(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                H = LennardJonesType::hessian(rij,r2,epsilon,sigma);
            }

            return H;
        }

    };

    using LennardJonesType1        = NonBondedHessian_<LennardJones_<BasicPotentials::LennardJones::Type1>>;
    using LennardJonesType2        = NonBondedHessian_<LennardJones_<BasicPotentials::LennardJones::Type2>>;
    using LennardJonesType3        = NonBondedHessian_<LennardJones_<BasicPotentials::LennardJones::Type3>>;

    using WCAType1                 = NonBondedHessian_<LennardJones_<BasicPotentials::WCA::Type1>>;
    using WCAType2                 = NonBondedHessian_<LennardJones_<BasicPotentials::WCA::Type2>>;
    using WCAType3                 = NonBondedHessian_<LennardJones_<BasicPotentials::WCA::Type3>>;

    using GeneralLennardJonesType1 = NonBondedHessian_<LennardJones_<BasicPotentials::GeneralLennardJones::Type1>>;
    using GeneralLennardJonesType2 = NonBondedHessian_<LennardJones_<BasicPotentials::GeneralLennardJones::Type2>>;
    using GeneralLennardJonesType3 = NonBondedHessian_<LennardJones_<BasicPotentials::GeneralLennardJones::Type3>>;

    using Steric6                  = NonBondedHessian_<LennardJones_<BasicPotentials::Steric::Steric6>>;
    using Steric12                 = NonBondedHessian_<LennardJones_<BasicPotentials::Steric::Steric12>>;

}}}}

#endif
