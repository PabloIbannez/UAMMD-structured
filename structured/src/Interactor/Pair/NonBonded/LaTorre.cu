#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicParameters/Pair/LaTorre.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    struct LaTorre_{

        using ParametersType        = typename BasicParameters::Pairs::LaTorre;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        struct ComputationalData{

            real4* pos;
            Box box;

            ParametersPairsIterator paramPairIterator;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> laTorreParam;
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

            computational.paramPairIterator = storage.laTorreParam->getPairIterator();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            ///////////////////////////////////////////////////////////

            storage.laTorreParam = std::make_shared<ParameterPairsHandler>(gd,pg,
                                                                      data);

            ///////////////////////////////////////////////////////////

            auto pairsParam = storage.laTorreParam->getPairParameters();

            real maxSigma = 0.0;
            real maxW     = 0.0;
            for(auto p : pairsParam){
                maxSigma=std::max(maxSigma,p.second.sigma);
                maxW    =std::max(maxW,p.second.w);
            }

            storage.cutOff = real(1.1224620)*maxSigma + maxW;

            System::log<System::MESSAGE>("[laTorre] cutOff: %f" ,storage.cutOff);

            return storage;
        }

        static inline __device__ real energy(int index_i, int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.paramPairIterator(index_i,index_j).epsilon;
            const real sigma   = computational.paramPairIterator(index_i,index_j).sigma;
            const real w       = computational.paramPairIterator(index_i,index_j).w;

            const real r2 = dot(rij, rij);
            const real r  = sqrt(r2);

            const real sgmFactor  = sigma*real(1.1224620);
            const real sgmFactorW = sgmFactor + w;

            real e;
            if        (r <= sgmFactor){
                e = -epsilon;
            } else if (r <= sgmFactorW){
                real x = real(M_PI)/(real(2.0)*w);
                     x = x*(r - sgmFactor);
                real cos2x = cos(x);
                     cos2x = cos2x*cos2x;
                     cos2x = std::min(cos2x,real(1.0));

                e = -epsilon*cos2x;
            } else { // r > sgmFactorW
                e = real(0.0);
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
            const real w       = computational.paramPairIterator(index_i,index_j).w;

            const real r2 = dot(rij, rij);
            const real r  = sqrt(r2);

            const real sgmFactor  = sigma*real(1.1224620);
            const real sgmFactorW = sgmFactor + w;

            real3 f;
            if        (r <= sgmFactor){
                f = make_real3(0.0);
            } else if (r <= sgmFactorW){
                const real x = real(M_PI)/(real(2.0)*w);
                const real y = x*(r - sgmFactor);
                      real sinycosy = sin(y)*cos(y);
                           sinycosy = std::min(sinycosy,real(1.0));
                           sinycosy = std::max(sinycosy,real(-1.0));
                const real factor = -real(2.0)*epsilon*x*(sinycosy);
                f = -factor*rij/r;
            } else { // r > sgmFactorW
                f = make_real3(0.0);
            }

            return f;
        }

    };

    using LaTorre  = NonBonded_<LaTorre_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,LaTorre,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::LaTorre>
)

