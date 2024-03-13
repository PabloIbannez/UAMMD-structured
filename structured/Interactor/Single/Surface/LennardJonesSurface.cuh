#ifndef __LENNARD_JONES_SURFACE_POT__
#define __LENNARD_JONES_SURFACE_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

    template<class LennardJonesSurfaceType>
    struct LennardJonesSurface_{

        using ParametersSingleType   = typename BasicParameters::Single::LennardJones;
        using ParameterSingleHandler = typename structured::SingleParameterHandler<ParametersSingleType>;

        using ParametersSingleIterator = typename ParameterSingleHandler::SingleIterator;

        struct ComputationalData{
            real4* pos;
            ParametersSingleIterator paramSingleIterator;
            real surfacePosition;
        };

        //Potential parameters
        struct StorageData{
            std::shared_ptr<ParameterSingleHandler> ljParam;
            real surfacePosition;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>                   gd,
                                                                 std::shared_ptr<ParticleGroup>  pg,
                                                                 DataEntry& data){

            StorageData storage;

            storage.ljParam = std::make_shared<ParameterSingleHandler>(gd,pg,
                                                                     data);

            storage.surfacePosition = data.getParameter<real>("surfacePosition",0);

            return storage;
        }


        //Computational data getter
        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>      gd,
                                                                 std::shared_ptr<ParticleGroup> pg,
                                                                 const StorageData&  storage){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.paramSingleIterator = storage.ljParam->getSingleIterator();

            computational.surfacePosition = storage.surfacePosition;

            return computational;
        }

        //Storage data reader


        static inline __device__  real3 force(const int& index_i,
                                              ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

            real3 p = make_real3(posi);
            real3 f = LennardJonesSurfaceType::force(p,computational.surfacePosition,
                                                     epsilon,sigma);

            return f;

        }

        static inline __device__  real energy(const int& index_i,
                                              ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

            real3 p = make_real3(posi);
            real  e = LennardJonesSurfaceType::energy(p,computational.surfacePosition,
                                                      epsilon,sigma);

            return e;

        }

      static inline __device__  real hessian(const int& index_i,
					     ComputationalData computational){

            real4 posi = computational.pos[index_i];

            const real epsilon = computational.paramSingleIterator(index_i).epsilon;
            const real sigma   = computational.paramSingleIterator(index_i).sigma;

            real3 p = make_real3(posi);
            real  e = LennardJonesSurfaceType::hessian(p,computational.surfacePosition,
						       epsilon,sigma);

            return e;

        }

    };

    using SurfaceLennardJonesType1        = Surface_<LennardJonesSurface_<BasicPotentials::Surface::LennardJones::Type1>>;
    using SurfaceLennardJonesType2        = Surface_<LennardJonesSurface_<BasicPotentials::Surface::LennardJones::Type2>>;

    using SurfaceWCAType1                 = Surface_<LennardJonesSurface_<BasicPotentials::Surface::WCA::Type1>>;
    using SurfaceWCAType2                 = Surface_<LennardJonesSurface_<BasicPotentials::Surface::WCA::Type2>>;

    using SurfaceGeneralLennardJonesType1 = Surface_<LennardJonesSurface_<BasicPotentials::Surface::GeneralLennardJones::Type1>>;
    using SurfaceGeneralLennardJonesType2 = Surface_<LennardJonesSurface_<BasicPotentials::Surface::GeneralLennardJones::Type2>>;

    using SurfaceAnchorage                = Surface_<LennardJonesSurface_<BasicPotentials::Surface::HarmonicWell>>;

}}}}

#endif
