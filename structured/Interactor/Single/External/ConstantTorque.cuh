#ifndef __EXTERNAL_CONSTANT_TORQUE_POT__
#define __EXTERNAL_CONSTANT_TORQUE_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct ConstantTorque_{

        //Computational data
        struct ComputationalData{
            real3 constantTorque;
        };

        struct StorageData{
            real3 constantTorque;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp){

            ComputationalData computational;
            computational.constantTorque = storage.constantTorque;
            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.constantTorque = data.getParameter<real3>("constantTorque");

            System::log<System::MESSAGE>("[ConstantTorque] Constant torque: %f %f %f",
                                         storage.constantTorque.x,
                                         storage.constantTorque.y,
                                         storage.constantTorque.z);

            return storage;
        }

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
            return real(0.0);
        }


      static inline __device__ ForceTorque forceTorque(int index_i,const ComputationalData& computational){
	ForceTorque forceTorque;
	forceTorque.force  = real4();
	forceTorque.torque = make_real4(computational.constantTorque, real(0.0));
	return forceTorque;
        }

    };

    using ConstantTorque   = ExternalTorque_<ConstantTorque_>;

}}}}

#endif
