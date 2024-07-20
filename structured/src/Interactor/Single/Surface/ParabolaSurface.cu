#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/Surface/Surface.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Surface.cuh"

#include "Interactor/BasicParameters/Single/Epsilon.cuh"

#include "Utils/ParameterHandler/SingleParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

    struct ParabolaSurface_{

        using ParametersSingleType   = typename BasicParameters::Single::Epsilon;
        using ParameterSingleHandler = typename structured::SingleParameterHandler<ParametersSingleType>;

        using ParametersSingleIterator = typename ParameterSingleHandler::SingleIterator;

        struct ComputationalData{
            real4* pos;
            ParametersSingleIterator paramSingleIterator;
            real surfacePosition;
        };

        //Potential parameters
        struct StorageData{
            std::shared_ptr<ParameterSingleHandler> epsilonParam;
            real surfacePosition;
        };

        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>      gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.paramSingleIterator = storage.epsilonParam->getSingleIterator();

            computational.surfacePosition = storage.surfacePosition;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                                 std::shared_ptr<ParticleGroup>        pg,
                                                                 DataEntry& data){

            StorageData storage;

            storage.epsilonParam = std::make_shared<ParameterSingleHandler>(gd,pg,
                                                                            data);

            storage.surfacePosition = data.getParameter<real>("surfacePosition",0);

            return storage;
        }


        static inline __device__  real3 force(const int& index_i,
                                              ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;

            real3 p = make_real3(posi);
            real3 f = BasicPotentials::Surface::Parabola::force(p,computational.surfacePosition,
                                                                  epsilon);

            return f;

        }

        static inline __device__  real energy(const int& index_i,
                                              ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;

            real3 p = make_real3(posi);
            real  e = BasicPotentials::Surface::Parabola::energy(p,computational.surfacePosition,
                                                                 epsilon);

            return e;

        }

    };

    using ParabolaSurface = Surface_<ParabolaSurface_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    Surface,ParabolaSurface,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::Surface::ParabolaSurface>
)

