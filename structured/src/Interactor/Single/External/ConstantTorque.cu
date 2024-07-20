#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/External/External.cuh"
#include "Interactor/InteractorFactory.cuh"

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
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

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

REGISTER_SINGLE_INTERACTOR(
    External,ConstantTorque,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::ConstantTorque>
)
