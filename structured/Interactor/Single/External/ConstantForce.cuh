#ifndef __EXTERNAL_CONSTANT_FORCE_POT__
#define __EXTERNAL_CONSTANT_FORCE_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct ConstantForce_{

        //Computational data
        struct ComputationalData{
            real3 constantForce;
        };

        struct StorageData{
            real3 constantForce;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){

            ComputationalData computational;
            computational.constantForce = storage.constantForce;
            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.constantForce = data.getParameter<real3>("constantForce");

            System::log<System::MESSAGE>("[ConstantForce] Constant force: %f %f %f",
                                         storage.constantForce.x,
                                         storage.constantForce.y,
                                         storage.constantForce.z);

            return storage;
        }

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
            return real(0.0);
        }


        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
            return computational.constantForce;
        }

    };

    using ConstantForce   = External_<ConstantForce_>;

}}}}

#endif
