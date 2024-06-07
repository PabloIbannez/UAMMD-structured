#ifndef __LINKER_POT__
#define __LINKER_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace SurfacePatches{

    struct Linker_{

        using ParametersSingleType   = typename BasicParameters::Single::LennardJones;
        using ParameterSingleHandler = typename structured::SingleParameterHandler<ParametersSingleType>;

        using ParametersSingleIterator = typename ParameterSingleHandler::SingleIterator;

        //Computational data
        struct ComputationalData{

            real4* patchesPos;
            real4* patchesVector;

            ParametersSingleIterator paramSingleIterator;

            real surfacePosition;
        };

        //Potential parameters
        struct StorageData{
            std::shared_ptr<ParameterSingleHandler> ljParam;
            real surfacePosition;
        };

        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                        std::shared_ptr<ParticleGroup> pg,
                                                        std::shared_ptr<GlobalData>    patchesGd,
                                                        std::shared_ptr<ParticleGroup> patchesPg,
                                                        DataEntry& data){

            StorageData storage;

            storage.ljParam = std::make_shared<ParameterSingleHandler>(patchesGd, patchesPg,
                                                                     data);

            storage.surfacePosition = data.getParameter<real>("surfacePosition",0);

            return storage;
        }

        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                      std::shared_ptr<ParticleGroup> pg,
                                                      std::shared_ptr<GlobalData>    patchesGd,
                                                      std::shared_ptr<ParticleGroup> patchesPg,
                                                      const StorageData& storage,
                                                      const Computables& comp,
                                                      const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();

            computational.patchesPos    = patchesPd->getPos(access::location::gpu, access::mode::read).raw();
            computational.patchesVector = patchesPd->getPatchVector(access::location::gpu, access::mode::read).raw();

            computational.paramSingleIterator = storage.ljParam->getSingleIterator();

            computational.surfacePosition = storage.surfacePosition;

            return computational;
        }

        static inline __device__ real energy(const int& index_i,
                                             const ComputationalData computational){

            const real4 patch_i_pos = computational.patchesPos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

            real e = real(0.0);

            e = BasicPotentials::Surface::GaussianWell::energy(make_real3(patch_i_pos),
                                                               computational.surfacePosition,
                                                               epsilon,sigma);

            return e;
        }


        static inline __device__ ForceTorque forceTorque(const int& index_i,
                                                         const ComputationalData computational){

            const real4 patch_i_pos = computational.patchesPos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

            ForceTorque frcTrq;

            real3 f  = make_real3(0.0);
            real3 t  = make_real3(0.0);

            f = BasicPotentials::Surface::GaussianWell::force(make_real3(patch_i_pos),
                                                              computational.surfacePosition,
                                                              epsilon,sigma);

            t = cross(make_real3(computational.patchesVector[index_i]),make_real3(f));

            frcTrq.force  = make_real4(f,0.0);
            frcTrq.torque = make_real4(t,0.0);

            return frcTrq;
        }

    };

    using Linker = SurfacePatchyParticles_<Linker_>;

}}}}

#endif
