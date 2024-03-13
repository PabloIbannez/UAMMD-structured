#ifndef __LENNARD_JONES_SOFT_CORE_POTENTIAL__
#define __LENNARD_JONES_SOFT_CORE_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template<class LennardJonesType>
    struct LennardJonesSoftCore_{

        using ParametersType        = typename BasicParameters::Pairs::LennardJones;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        struct ComputationalData{

            real4* pos;
            Box box;

            ParametersPairsIterator paramPairIterator;

            real cutOffFactor;

            real lambda;
            real alpha;
            int  n;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> ljParam;

            real cutOffFactor;
            real cutOff;

            real alpha;
            int  n;
        };

        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box    = gd->getEnsemble()->getBox();
            computational.lambda = gd->getEnsemble()->getLambda();

            computational.paramPairIterator = storage.ljParam->getPairIterator();

            computational.cutOffFactor = storage.cutOffFactor;
            computational.alpha        = storage.alpha;
            computational.n            = storage.n;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;
            storage.cutOffFactor  = data.getParameter<real>("cutOffFactor");
            storage.alpha         = data.getParameter<real>("alpha");
            storage.n             = data.getParameter<int>("n", 2);

            System::log<System::MESSAGE>("[LenardJonesSoftCore] cutOffFactor: %f", storage.cutOffFactor);
            System::log<System::MESSAGE>("[LenardJonesSoftCore] alpha: %f", storage.alpha);
            System::log<System::MESSAGE>("[LenardJonesSoftCore] n: %i", storage.n);

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

            const real lambda  = computational.lambda;
            if(lambda==real(0.0)){
                return real(0.0);
            }

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            const real r2 = dot(rij, rij);

            real e = real(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                e = LennardJonesType::energy(rij,r2,epsilon,sigma,lambda,alpha,n);
            }

            return e;

        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             const ComputationalData& computational){

            const real lambda  = computational.lambda;
            if(lambda==real(0.0)){
                return make_real3(0.0);
            }

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            const real r2 = dot(rij, rij);

            real3 f = make_real3(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                f = LennardJonesType::force(rij,r2,epsilon,sigma,lambda,alpha,n);
            }

            return f;
        }

        static inline __device__ real lambdaDerivative(int index_i, int index_j,
                                                       const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;
            const real lambda  = computational.lambda;
            const real alpha   = computational.alpha;
            const int  n       = computational.n;

            const real r2 = dot(rij, rij);

            real ld = real(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                ld = LennardJonesType::lambdaDerivative(rij,r2,epsilon,sigma,lambda,alpha,n);
            }

            return ld;

        }

    };

    using LennardJonesSoftCoreType1 = NonBondedLambda_<LennardJonesSoftCore_<BasicPotentials::LennardJones::SoftCoreType1>>;
    using LennardJonesSoftCoreType2 = NonBondedLambda_<LennardJonesSoftCore_<BasicPotentials::LennardJones::SoftCoreType2>>;

    using Steric6SoftCore           = NonBondedLambda_<LennardJonesSoftCore_<BasicPotentials::Steric::Steric6SoftCore>>;
    using Steric12SoftCore          = NonBondedLambda_<LennardJonesSoftCore_<BasicPotentials::Steric::Steric12SoftCore>>;

}}}}

#endif
