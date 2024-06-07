#ifndef __EXTERNAL_AC_MAGNETICFIELD_POT__
#define __EXTERNAL_AC_MAGNETICFIELD_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

  struct ConstantMagneticField_{

    //Computational data
    struct ComputationalData{
      real4* dir;
      real4* magnetization;

      real3  magneticField;
    };

    struct StorageData{
      real  b0;
      real3 direction;
    };

    static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
							   std::shared_ptr<ParticleGroup> pg,
                                                           const StorageData&  storage,
                                                           const Computables& comp,
                                                           const cudaStream_t& st){

      ComputationalData computational;
      std::shared_ptr<ParticleData> pd = pg->getParticleData();

      computational.magnetization = pd->getMagnetization(access::location::gpu,access::mode::read).raw();
      computational.dir = pd->getDir(access::location::gpu,access::mode::read).raw();

      real b0 = storage.b0;
      real3 direction = storage.direction;

      computational.magneticField = b0*direction;

      return computational;
    }

    static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
						 std::shared_ptr<ParticleGroup> pg,
						 DataEntry& data){

      StorageData storage;

      storage.b0        = data.getParameter<real>("b0");
      storage.direction = data.getParameter<real3>("direction");

      System::log<System::MESSAGE>("[ConstantMagneticField] Amplitude, b0 = %f", storage.b0);
      System::log<System::MESSAGE>("[ConstantMagneticField] Direction, direction = (%f, %f, %f)",
                                   storage.direction.x, storage.direction.y, storage.direction.z);

      return storage;
    }

    static inline __device__ real energy(int index_i, const ComputationalData& computational){

      const real4 diri = computational.dir[index_i];
      const real4 m_and_M = computational.magnetization[index_i];
      real3 magneticMoment = m_and_M.w*rotateVector(diri, make_real3(m_and_M));
      real3 magneticField = computational.magneticField;
      real e = dot(magneticMoment, magneticField);
      return e;
    }

    static inline __device__ real4 magneticField(const int index_i,const ComputationalData& computational){
	  	return make_real4(computational.magneticField, 0);
	  }

    static inline __device__ ForceTorque forceTorque(const int index_i,const ComputationalData& computational){
      ForceTorque forceTorque;

      forceTorque.force  = make_real4(0.0);

      const real4 diri    = computational.dir[index_i];
      const real4 m_and_M = computational.magnetization[index_i];

      real3 magneticMoment = m_and_M.w*rotateVector(diri, make_real3(m_and_M));

      forceTorque.torque = make_real4(cross(magneticMoment, make_real3(magneticField(index_i, computational))), 0);

	  	return forceTorque;
	  }

  };

  //B(t) = b0*sin(w*t+phase); w = 2*pi*f
  struct ACMagneticField_ : public ConstantMagneticField_{

    using ConstantMagneticField_::ComputationalData;

    struct StorageData: public ConstantMagneticField_::StorageData{
      real frequency;
      real phase;
    };

    static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
	                                                   std::shared_ptr<ParticleGroup> pg,
                                                           const StorageData&  storage,
                                                           const Computables& comp,
                                                           const cudaStream_t& st){

      ComputationalData computational =
      ConstantMagneticField_::getComputationalData(gd, pg, storage, comp, st);

      real b0 = storage.b0;
      real w  = storage.frequency*2*M_PI;
      real phase = storage.phase;
      real time = (gd->getFundamental()->getCurrentStep())*(gd->getFundamental()->getTimeStep());
      real3 direction = storage.direction;

      computational.magneticField = b0*sin(w*time+phase)*direction;

      return computational;
    }

    static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                  						 std::shared_ptr<ParticleGroup> pg,
                                  						 DataEntry& data){

      StorageData storage;
      static_cast<ConstantMagneticField_::StorageData&>(storage) = ConstantMagneticField_::getStorageData(gd, pg, data);

      storage.frequency = data.getParameter<real>("frequency");
      storage.phase     = data.getParameter<real>("phase", 0);

      System::log<System::MESSAGE>("[ACMagneticField] Frequency = %f", storage.frequency);
      System::log<System::MESSAGE>("[ACMagneticField] Phase = %f", storage.phase);
      System::log<System::MESSAGE>("[ACMagneticField] Phase = %f %f %f",
			                             storage.direction.x, storage.direction.y, storage.direction.z);

      return storage;
    }

    static inline __device__ real energy(int index_i, const ComputationalData& computational){
      return ConstantMagneticField_::energy(index_i, computational);
    }

    static inline __device__ real4 magneticField(const int index_i,const ComputationalData& computational){
      return ConstantMagneticField_::magneticField(index_i, computational);
    }

    static inline __device__ ForceTorque forceTorque(const int index_i,const ComputationalData& computational){
      return ConstantMagneticField_::forceTorque(index_i, computational);
    }

  };

  using ConstantMagneticField = ExternalForceTorqueMagneticField_<ConstantMagneticField_>;
  using ACMagneticField       = ExternalForceTorqueMagneticField_<ACMagneticField_>;

}}}}

#endif
