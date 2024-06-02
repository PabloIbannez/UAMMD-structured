#ifndef __EXTERNAL_CONSTANT_MAGNETICFIELD_POT__
#define __EXTERNAL_CONSTANT_MAGNETICFIELD_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct ConstantMagneticField_{

        //Computational data
        struct ComputationalData{
      	  real4* dir;
      	  real4* magnetization;
      	  real3 constantMagneticField;
        };

        struct StorageData{
            real3 constantMagneticField;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;
	    std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.magnetization = pd->getMagnetization(access::location::gpu,
							       access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu,
					   access::mode::read).raw();
            computational.constantMagneticField = storage.constantMagneticField;
            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.constantMagneticField = data.getParameter<real3>("constantMagneticField");

            System::log<System::MESSAGE>("[ConstantMagneticField] Constant magneticField: %f %f %f",
                                         storage.constantMagneticField.x,
                                         storage.constantMagneticField.y,
                                         storage.constantMagneticField.z);

            return storage;
        }

//      static inline __device__ ForceTorqueMagneticField forceTorqueMagneticField(int index_i, const ComputationalData& computational){
//	ForceTorqueMagneticField f_t_mf;
//	f_t_mf.force = real4();
//	const real4 diri    = computational.dir[index_i];
//	const real4 m_and_M = computational.magnetization[index_i];
//
//	real3 magneticMoment = m_and_M.w*rotateVector(diri, make_real3(m_and_M));
//	real3 magneticField  = computational.constantMagneticField;
//
//	f_t_mf.magneticField = make_real4(magneticField, 0);
//	f_t_mf.torque        = make_real4(cross(magneticMoment, magneticField), 0);
//	return f_t_mf;
//      }

        static inline __device__ real energy(int index_i, const ComputationalData& computational){

	  const real4 diri = computational.dir[index_i];
	  const real4 m_and_M = computational.magnetization[index_i];
	  real3 magneticMoment = m_and_M.w*rotateVector(diri, make_real3(m_and_M));
	  real3 magneticField = computational.constantMagneticField;
	  real e = dot(magneticMoment, magneticField);
	  return e;
        }

      static inline __device__ ForceTorque forceTorque(const int index_i,const ComputationalData& computational){
        ForceTorque forceTorque;

        forceTorque.force  = make_real4(0.0);
        forceTorque.torque = make_real4(0.0);

  	  	return forceTorque;
  	  }

  	  static inline __device__ real4 magneticField(const int index_i,const ComputationalData& computational){

  	  	real4 magneticField = make_real4(0.0);

  	  	return magneticField;
  	  }

    };


    using ConstantMagneticField = ExternalTorqueMagneticField_<ConstantMagneticField_>;

}}}}

#endif
