#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Surface.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

    struct Absorbed_{

        struct ComputationalData{

            Box box;

            real4* pos;
            int*   id;

            const int* id2index;
            real4* absorbedPos;
            bool*  isAbsorbed;

            real K;
            real heightThreshold;
        };

        //Potential parameters
        struct StorageData{

            thrust::device_vector<real4> absorbedPos;
            thrust::device_vector<bool>  isAbsorbed;

            real K;
            real heightThreshold;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K               = data.getParameter<real>("K");
            storage.heightThreshold = data.getParameter<real>("heightThreshold");

            int N = pg->getParticleData()->getNumParticles();

            storage.absorbedPos.resize(N);
            storage.isAbsorbed.resize(N);

            return storage;
        }


        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            static bool computedAbsorbed = false;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            if(not computedAbsorbed){

                System::log<System::MESSAGE>("[Absorbed] Computing absorbed particles.");

                thrust::host_vector<real4> absorbedPos_h;
                thrust::host_vector<bool>  isAbsorbed_h;

                int N = pg->getParticleData()->getNumParticles();

                absorbedPos_h.resize(N);
                isAbsorbed_h.resize(N);

                auto ids = pd->getId(access::location::cpu,access::mode::read);
                auto pos = pd->getPos(access::location::cpu,access::mode::read);

                int absorbedCount = 0;
                fori(0,N){
                    int   id = ids[i];
                    real4 p  = pos[i];

                    if(p.z < storage.heightThreshold){
                        absorbedPos_h[id] = p;
                        isAbsorbed_h[id]  = true;
                        absorbedCount++;
                    }
                    else{
                        isAbsorbed_h[id]  = false;
                    }

                }

                System::log<System::MESSAGE>("[Absorbed] %d particles absorbed.",absorbedCount);

                storage.absorbedPos = absorbedPos_h;
                storage.isAbsorbed  = isAbsorbed_h;

                computedAbsorbed = true;

                System::log<System::MESSAGE>("[Absorbed] Absorbed particles computed.");
            }

            ComputationalData computational;
            {
                computational.box = gd->getEnsemble()->getBox();

                computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
                computational.id  = pd->getId(access::location::gpu,  access::mode::read).raw();

                computational.id2index = pd->getIdOrderedIndices(access::location::gpu);
                computational.absorbedPos = thrust::raw_pointer_cast(storage.absorbedPos.data());
                computational.isAbsorbed  = thrust::raw_pointer_cast(storage.isAbsorbed.data());

                computational.K               = storage.K;
                computational.heightThreshold = storage.heightThreshold;
            }

            return computational;
        }

        //Storage data reader
        static inline __device__  real3 force(const int& index_i,
                                              ComputationalData computational){

            real3 f = make_real3(0.0);

            const int id_i = computational.id[index_i];

            if(computational.isAbsorbed[id_i]){

                real3 posi = make_real3(computational.pos[index_i]);
                real3 posj = make_real3(computational.absorbedPos[id_i]);

                real3 rij = computational.box.apply_pbc(posj - posi);

                real3 K   = make_real3(computational.K);
                real3 r0  = make_real3(0.0);

                f = BasicPotentials::HarmonicAnisotropic::force(rij,K,r0);
            }

            return f;

        }

        static inline __device__  real energy(const int& index_i,
                                              ComputationalData computational){

            real e = real(0.0);

            const int id_i = computational.id[index_i];

            if(computational.isAbsorbed[id_i]){
                real3 posi = make_real3(computational.pos[index_i]);
                real3 posj = make_real3(computational.absorbedPos[id_i]);

                real3 rij = computational.box.apply_pbc(posj - posi);

                real3 K   = make_real3(computational.K);
                real3 r0  = make_real3(0.0);

                e = BasicPotentials::HarmonicAnisotropic::energy(rij,K,r0);
            }

            return e;

        }

    };

    using Absorbed = Surface_<Absorbed_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    Surface,Absorbed,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::Surface::Absorbed>
)
