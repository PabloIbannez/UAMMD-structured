#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{


    struct Clashed_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real*  radius;

            Box box;

            real lambda;
            real gamma;
        };

        //Potential parameters
        struct StorageData{

            real lambda;
            real gamma;

            real cutOff;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                std::shared_ptr<ParticleGroup> pg,
                const StorageData& storage,
                const Computables& comp,
                const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.lambda = storage.lambda;
            computational.gamma  = storage.gamma;

            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                std::shared_ptr<ParticleGroup> pg,
                DataEntry& data){

            StorageData storage;

            storage.lambda  = data.getParameter<real>("lambda");
            storage.gamma   = data.getParameter<real>("gamma");

            //////////////////////////////////////

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            auto radius     = pd->getRadius(access::location::cpu, access::mode::readwrite);
            auto groupIndex = pg->getIndexIterator(access::location::cpu);

            real maxRadius = radius[0];
            fori(0,pg->getNumberParticles()){
                int  index = groupIndex[i];
                maxRadius=std::max(maxRadius,radius[index]);
            }

            //////////////////////////////////////

            storage.cutOff = real(2.0)*maxRadius*storage.gamma;

            System::log<System::MESSAGE>("[Clashed] cutOff: %f" ,storage.cutOff);

            return storage;
        }

        static inline __device__ real energy(int index_i, int index_j,
                const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real d  = computational.gamma*(computational.radius[index_i]+computational.radius[index_j]);
            real d2 = d*d;

            real e = max(real(0.0),d2-r2);
            e = e*e;

            return computational.lambda*e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                const ComputationalData& computational){

            const real4 posi =  computational.pos[index_i];
            const real4 posj =  computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real d  = computational.gamma*(computational.radius[index_i]+computational.radius[index_j]);
            real d2 = d*d;

            real fmod = 0;

            if(r2>0 and d2>=r2){
                fmod = -real(4.0)*(d2-r2); //fmod*rij/r=-real(4.0)*r*(d2-r2)*rij/r
            }

            return computational.lambda*fmod*rij;
        }

        static inline __device__ tensor3 hessian(int index_i, int index_j,
                const ComputationalData& computational){

            const real4 posi =  computational.pos[index_i];
            const real4 posj =  computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);
            const real r    = sqrt(r2);

            real d  = computational.gamma*(computational.radius[index_i]+computational.radius[index_j]);
            real d2 = d*d;

            tensor3 H = tensor3();

            if(r2>0 and d2>=r2){
                const real invr2 = real(1.0)/r2;
                const real invr  = real(1.0)/r;

                real energyDerivative       = -real(4.0)*(d2-r2)*r;
                real energySecondDerivative = -real(4.0)*(d2-real(3.0)*r2);

                H = computeHessianRadialPotential(rij, invr, invr2,
                        energyDerivative,
                        energySecondDerivative);
            }

            return computational.lambda*H;
        }

    };

    using Clashed = NonBonded_<Clashed_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,Clashed,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::Clashed>
)
