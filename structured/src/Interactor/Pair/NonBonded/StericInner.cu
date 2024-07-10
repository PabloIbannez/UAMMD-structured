#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Steric.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template<class StericType>
    struct StericInner_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real* innerRadius;

            Box   box;

            real epsilon;
            real cutOffFactor;
        };

        //Potential parameters
        struct StorageData{

            real epsilon;

            real cutOffFactor;
            real cutOff;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.epsilon      = data.getParameter<real>("epsilon");
            storage.cutOffFactor = data.getParameter<real>("cutOffFactor");

            ///////////////////////////////////////////////////////////

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            auto innerRadiusIterator = pg->getPropertyIterator(pd->getInnerRadius(access::location::cpu, access::mode::read).begin(),
                                                               access::location::cpu);

            int N = pg->getNumberParticles();

            std::vector<real> sortedRadius(innerRadiusIterator, innerRadiusIterator + N);

            //Sort in descending order
            std::sort(sortedRadius.begin(), sortedRadius.end(), std::greater<real>());
            storage.cutOff = (sortedRadius[0] + sortedRadius[1]) * storage.cutOffFactor;

            System::log<System::MESSAGE>("[StericInner] cutOff: %f" ,storage.cutOff);

            return storage;
        }

        ///////////////////////////

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos         = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.innerRadius = pd->getInnerRadius(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.epsilon      = storage.epsilon;
            computational.cutOffFactor = storage.cutOffFactor;

            return computational;
        }

        static inline __device__ real energy(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.epsilon;
            const real sigma   = computational.innerRadius[index_i] + computational.innerRadius[index_j];

            const real r2 = dot(rij, rij);

            real e = real(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                e = StericType::energy(rij,r2,epsilon,sigma);
            }

            return e;

        }

        static inline __device__ real3 force(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));

            const real epsilon = computational.epsilon;
            const real sigma   = computational.innerRadius[index_i] + computational.innerRadius[index_j];

            const real r2 = dot(rij, rij);

            real3 f = make_real3(0.0);

            real cutOff2=sigma*computational.cutOffFactor;
                 cutOff2=cutOff2*cutOff2;
            if(r2<=cutOff2){
                f = StericType::force(rij,r2,epsilon,sigma);
            }

            return f;
        }
    };

    using StericInner6  = NonBonded_<StericInner_<BasicPotentials::Steric::Steric6>>;
    using StericInner12 = NonBonded_<StericInner_<BasicPotentials::Steric::Steric12>>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,StericInner6,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::StericInner6>
)
