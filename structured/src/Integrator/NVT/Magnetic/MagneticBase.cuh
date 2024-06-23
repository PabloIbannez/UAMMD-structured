#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

#include "Landau_Lifshitz_Gilbert.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace Magnetic{

  template<class T>  using cached_vector = uninitialized_cached_vector<T>;
  class MagneticBase : public Integrator{

  private:

    cudaStream_t stream;
    std::shared_ptr<GlobalData> gd;

    //////////////////////////////

    real dt;
    real kBT;
    real damping;
    real msat;
    real gyroRatio;
    //////////////////////////////
    //Write here any defined magnetic integration algorithm
    const std::vector<std::string> validMagneticIntegrationAlgorithms = {"none", "LLG_Euler", "LLG_Heun"};
    std::string magneticIntegrationAlgorithm;

  public:

    struct Parameters{
      cudaStream_t stream;
      real dt;
      real kBT;
      real damping;
      real msat;
      real gyroRatio;
      std::string magneticIntegrationAlgorithm;
    };

    MagneticBase(std::shared_ptr<GlobalData>    gd,
		 std::shared_ptr<ParticleGroup> pg,
		 Parameters             parameters,
		 std::string name):Integrator(pg,name),gd(gd){

      System::log<System::MESSAGE>("[MagneticBase] Created MagneticBase integrator \"%s\"",name.c_str());

      stream = parameters.stream;
      dt = parameters.dt;
      kBT = parameters.kBT;
      damping = parameters.damping;
      msat = parameters.msat;
      gyroRatio = parameters.gyroRatio;
      if (checkValidMagneticIntegrationAlgorithm(parameters.magneticIntegrationAlgorithm)){
	magneticIntegrationAlgorithm = parameters.magneticIntegrationAlgorithm;
	System::log<System::MESSAGE>("[MagneticBase] Using the \"%s\" algorithm to integrate the internal dynamics of the magnetization", this->magneticIntegrationAlgorithm.c_str());
      } else {
	if (parameters.magneticIntegrationAlgorithm == "NotSelected")
	  System::log<System::CRITICAL>("[MagneticBase] Please select a valid integration algorithm");
	else
	  System::log<System::CRITICAL>("[MagneticBase] The selected magnetic integration algorithm is not valid");
      }
      if (parameters.magneticIntegrationAlgorithm != "none"){
	if (damping<0)
	  System::log<System::CRITICAL>("[MagneticBase] Please, set the Gilbert damping parameter");
	if (gyroRatio<0)
	  System::log<System::CRITICAL>("[MagneticBase] Please, set the gyromagnetic ratio");
      }
    }

    MagneticBase(std::shared_ptr<GlobalData>    gd,
		 std::shared_ptr<ParticleGroup> pg,
		 DataEntry&                   data,
		 std::string name):Integrator(pg,name),gd(gd){

      Parameters parameters;
      parameters.stream = this->stream;
      parameters.dt     = this->dt;
      parameters.kBT		= this->kBT;
      parameters.damping = data.getParameter<real>("damping");
      parameters.msat = data.getParameter<real>("msat");
      parameters.gyroRatio = data.getParameter<real>("gyroRatio");
      parameters.magneticIntegrationAlgorithm = data.getParameter<std::string>("magneticIntegrationAlgorithm","NotSelected");
      MagneticBase(gd, pg, parameters, name);
    }

      void resetMagneticField(){

          auto magneticField = pd->getMagneticField(access::location::gpu, access::mode::readwrite);
          thrust::fill(thrust::cuda::par.on(stream), magneticField.begin(), magneticField.end(), make_real4(0));

          //CudaSafeCall(cudaStreamSynchronize(stream));
          CudaSafeCall(cudaDeviceSynchronize());
          CudaCheckError();
      }

      void updateMagneticField() {
          for(auto forceComp: interactors) forceComp->sum({.force =false,
                                                           .energy=false,
                                                           .virial=false,
                                                           .stress=false,
                                                           .magneticField=true},stream);
          //CudaSafeCall(cudaStreamSynchronize(stream));
          CudaSafeCall(cudaDeviceSynchronize());
          CudaCheckError();
      }


	  void forwardTime() override {
	    //Warning: note simulation step is not updated !!!
	    updateMagneticField();
	    updateMagnetization();
	    resetMagneticField();
	  }

	  void updateMagnetization(){
	    if (magneticIntegrationAlgorithm == "LLG_Euler"){

	      LLG::Euler::updateMagnetization(pd, pg, gd, stream,
					      damping, msat, gyroRatio,
					      dt, kBT);

	    } else if (magneticIntegrationAlgorithm == "LLG_Heun"){
	      cached_vector<real3> initialMagnetization(pg->getNumberParticles());
	      LLG::Heun::updateHalfStep(pd, pg, gd, initialMagnetization,
					damping, msat, gyroRatio,
					dt, kBT, 0, stream);
	      resetMagneticField();
	      CudaCheckError();
	      real timeStep = gd->getFundamental()->getTimeStep();
	      int currentStep = gd->getFundamental()->getCurrentStep();
	      //I need to upload the time to t+1/2 = dt*(currentStep+0.5), but I can't add 0.5 to
	      //currentStep because it is integer. For that reason The solution is change dt to
	      //a value such that the time is t+1/2.
	      //I undo this change after computing the magnetic fields
	      //gd->setTimeStep(timeStep*(1+0.5/currentStep));
	      updateMagneticField();
	      //gd->setTimeStep(timeStep);
	      LLG::Heun::updateHalfStep(pd, pg, gd, initialMagnetization,
					damping, msat, gyroRatio,
					dt, kBT, 1, stream);
	    }
	  }

  private:

    bool checkValidMagneticIntegrationAlgorithm(const std::string algorithm_in){
      return std::find(validMagneticIntegrationAlgorithms.begin(),
		       validMagneticIntegrationAlgorithms.end(),
		       algorithm_in) != validMagneticIntegrationAlgorithms.end();
    }

  };

}}}}}

