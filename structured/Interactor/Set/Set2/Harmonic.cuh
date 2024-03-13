#ifndef __SET2_HARMONIC__
#define __SET2_HARMONIC__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set2{

    struct Harmonic_{

        struct ComputationalData{
            real4* pos;
            real* mass;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct SetParameters{
            real K;
            real r0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                                 std::shared_ptr<ExtendedParticleData> pd,
                                                                 const StorageData&  storage){
            ComputationalData computational;

            computational.pos  = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.mass = pd->getMass(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ExtendedParticleData> pd,
                                                   DataEntry& data){

            StorageData storage;
            return storage;
        }

        //Set parameters reader

        template<typename T>
        static __host__ SetParameters processSetParameters(std::shared_ptr<GlobalData> gd,
                                                           std::map<std::string,T>& bondParametersMap){

            SetParameters param;

            param.K  = bondParametersMap.at("K");
            param.r0 = bondParametersMap.at("r0");

            return param;
        }

        static inline __device__ real energy(const real3& centerOfMass1,
                                             const real3& centerOfMass2,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 r12 = computational.box.apply_pbc(centerOfMass2-centerOfMass1);
            const real  r2  = dot(r12,r12);

            return BasicPotentials::Harmonic::energy(r12,r2,setParam.K,setParam.r0);
        }

        static inline __device__ real4 force(const real3& centerOfMass1,
                                             const real3& centerOfMass2,
                                             const ComputationalData& computational,
                                             const SetParameters& setParam){

            const real3 r12 = computational.box.apply_pbc(centerOfMass2-centerOfMass1);
            const real  r2  = dot(r12,r12);

            return make_real4(BasicPotentials::Harmonic::force(r12,r2,setParam.K,setParam.r0));
        }

    };

    using HarmonicBondBetweenCentersOfMass  = CenterOfMass_<Harmonic_>;

}}}}

#endif
