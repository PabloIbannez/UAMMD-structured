#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"
#include "Interactor/BasicParameters.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template<class LennardJonesType>
    struct DLVO_{

        using ParametersType        = typename BasicParameters::Pairs::LennardJones;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

          real4* pos;
          real*  radius;
          real*  charge;

          Box box;

          ParametersPairsIterator paramPairIterator;

          real ELECOEF;

          real dielectricConstant;
          real debyeLength;

          real cutOffNPFactor;
          real cutOffDHFactor;
        };

        //Potential parameters
        struct StorageData{

          std::shared_ptr<ParameterPairsHandler> ljParam;

          real ELECOEF;

          real dielectricConstant;
          real debyeLength;

          real cutOffNPFactor;
          real cutOffDHFactor;

          real cutOff;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();
            computational.charge = pd->getCharge(access::location::gpu, access::mode::read).raw();

            computational.box    = gd->getEnsemble()->getBox();

            computational.paramPairIterator = storage.ljParam->getPairIterator();

            computational.ELECOEF = storage.ELECOEF;

            computational.dielectricConstant = storage.dielectricConstant;
            computational.debyeLength = storage.debyeLength;

            computational.cutOffNPFactor = storage.cutOffNPFactor;
            computational.cutOffDHFactor = storage.cutOffDHFactor;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.ELECOEF = gd->getUnits()->getElectricConversionFactor();

            storage.cutOffNPFactor  = data.getParameter<real>("cutOffNPFactor");
            storage.cutOffDHFactor  = data.getParameter<real>("cutOffDHFactor");

            storage.dielectricConstant = data.getParameter<real>("dielectricConstant");
            storage.debyeLength        = data.getParameter<real>("debyeLength");

            ///////////////////////////////////////////////////////////

            storage.ljParam = std::make_shared<ParameterPairsHandler>(gd, pg,
                                                                    data);

            ///////////////////////////////////////////////////////////

            auto pairsParam = storage.ljParam->getPairParameters();

            real maxSigma = 0.0;
            for(auto p : pairsParam){
                maxSigma=std::max(maxSigma,p.second.sigma);
            }

            real cutOffNP = maxSigma*storage.cutOffNPFactor;
            real cutOffDH = storage.debyeLength*storage.cutOffDHFactor;

            storage.cutOff = std::max(cutOffNP,cutOffDH);

            System::log<System::MESSAGE>("[DLVO] cutOffNP: %f" ,cutOffNP);
            System::log<System::MESSAGE>("[DLVO] cutOffDH: %f" ,cutOffDH);

            return storage;
        }

        static inline __device__ real energy(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real e = real(0.0);

            const real sigma = computational.paramPairIterator(index_i, index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real epsilon = computational.paramPairIterator(index_i, index_j).epsilon;

                e += LennardJonesType::energy(rij,r2,epsilon,sigma);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                e += BasicPotentials::DebyeHuckel::DebyeHuckelSpheres::energy(rij,r2,computational.ELECOEF,chgProduct,
                                                                               computational.radius[index_i],computational.radius[index_j],
                                                                               computational.dielectricConstant,
                                                                               computational.debyeLength);
            }

            return e;

        }

        static inline __device__ real3 force(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real3 f = make_real3(0.0);

            const real sigma = computational.paramPairIterator(index_i, index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real epsilon = computational.paramPairIterator(index_i, index_j).epsilon;

                f += LennardJonesType::force(rij,r2,epsilon,sigma);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                f += BasicPotentials::DebyeHuckel::DebyeHuckelSpheres::force(rij,r2,computational.ELECOEF,chgProduct,
                                                                             computational.radius[index_i],computational.radius[index_j],
                                                                             computational.dielectricConstant,computational.debyeLength);
            }

            return f;
        }

      static inline __device__ tensor3 hessian(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            tensor3 H = tensor3(0.0);

            const real sigma = computational.paramPairIterator(index_i, index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real epsilon = computational.paramPairIterator(index_i, index_j).epsilon;

                H += LennardJonesType::hessian(rij,r2,epsilon,sigma);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                H += BasicPotentials::DebyeHuckel::DebyeHuckelSpheres::hessian(rij,r2,computational.ELECOEF,chgProduct,
									       computational.radius[index_i],computational.radius[index_j],
									       computational.dielectricConstant,computational.debyeLength);
            }

            return H;
        }

    };

    using DLVOType1 = NonBondedHessian_<DLVO_<BasicPotentials::LennardJones::Type1>>;
    using DLVOType2 = NonBondedHessian_<DLVO_<BasicPotentials::LennardJones::Type2>>;
    using DLVOType3 = NonBondedHessian_<DLVO_<BasicPotentials::LennardJones::Type3>>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,DLVOType1,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::DLVOType1>
)

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,DLVOType2,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::DLVOType2>
)

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,DLVOType3,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::DLVOType3>
)
