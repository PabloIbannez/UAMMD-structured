#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/DebyeHuckel.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    struct DebyeHuckel_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real*  charge;

            Box    box;

            real ELECOEF;

            real dielectricConstant;
            real debyeLength;

            real cutOffFactor;
        };

        //Potential parameters
        struct StorageData{

            real ELECOEF;

            real dielectricConstant;
            real debyeLength;

            real cutOffFactor;
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
            computational.charge = pd->getCharge(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.ELECOEF = storage.ELECOEF;

            computational.dielectricConstant = storage.dielectricConstant;
            computational.debyeLength = storage.debyeLength;

            computational.cutOffFactor = storage.cutOffFactor;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.ELECOEF = gd->getUnits()->getElectricConversionFactor();

            storage.dielectricConstant = data.getParameter<real>("dielectricConstant");
            storage.debyeLength        = data.getParameter<real>("debyeLength");

            storage.cutOffFactor = data.getParameter<real>("cutOffFactor");
            storage.cutOff       = storage.cutOffFactor*storage.debyeLength;

            System::log<System::MESSAGE>("[DebyeHuckel] cutOff: %f" ,storage.cutOff);

            return storage;

        }


        static inline __device__ real energy(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real e = real(0.0);

            real cutOff2 = computational.debyeLength*computational.cutOffFactor;
                 cutOff2 = cutOff2*cutOff2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOff2 and chgProduct != real(0.0)){

                e+=BasicPotentials::DebyeHuckel::DebyeHuckel::energy(rij,r2,
                                                                     computational.ELECOEF,
                                                                     chgProduct,
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

            real3 f = make_real3(real(0.0));

            real cutOff2 = computational.debyeLength*computational.cutOffFactor;
                 cutOff2 = cutOff2*cutOff2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOff2 and chgProduct != real(0.0)){

                f+=BasicPotentials::DebyeHuckel::DebyeHuckel::force(rij,r2,
                                                                    computational.ELECOEF,
                                                                    chgProduct,
                                                                    computational.dielectricConstant,
                                                                    computational.debyeLength);
            }

            return f;
        }

      static inline __device__ tensor3 hessian(const int index_i,const int index_j,
					       const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            tensor3 H = tensor3(real(0.0));

            real cutOff2 = computational.debyeLength*computational.cutOffFactor;
                 cutOff2 = cutOff2*cutOff2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOff2 and chgProduct != real(0.0)){

                H+=BasicPotentials::DebyeHuckel::DebyeHuckel::hessian(rij,r2,
								      computational.ELECOEF,
								      chgProduct,
								      computational.dielectricConstant,
								      computational.debyeLength);
            }

            return H;
        }
    };

    using DebyeHuckel = NonBonded_<DebyeHuckel_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,DebyeHuckel,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::DebyeHuckel>
)
