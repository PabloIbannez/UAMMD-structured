#include "MagneticBase.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace Magnetic{

  class Fixed : public IntegratorBaseNVT{

  private:

    private:

    std::unique_ptr<MagneticBase> magnetic;
    bool firstStep = true;

  public:

    Fixed(std::shared_ptr<GlobalData>           gd,
	  std::shared_ptr<ParticleGroup>        pg,
	  DataEntry& data,
	  std::string name):IntegratorBaseNVT(gd,pg,data,name){

      typename MagneticBase::Parameters parameters;

      parameters.stream = this->stream;
      parameters.dt     = this->dt;
      parameters.kBT		= this->kBT;
      parameters.damping = data.getParameter<real>("damping");
      parameters.msat = data.getParameter<real>("msat");
      parameters.gyroRatio = data.getParameter<real>("gyroRatio");
      parameters.magneticIntegrationAlgorithm = data.getParameter<std::string>("magneticIntegrationAlgorithm","NotSelected");

      magnetic = std::make_unique<MagneticBase>(gd,pg,parameters,name);

    }

    void forwardTime() override {

      if(firstStep){
	//Load all interactors into magnetic
	for(auto& interactor : this->getInteractors()){
	  magnetic->addInteractor(interactor);
	}

	//Load all updatables into magnetic
	for(auto& updatable : this->getUpdatables()){
	  magnetic->addUpdatable(updatable);
	}

	firstStep = false;
      }

      magnetic->updateMagneticField();
      magnetic->updateMagnetization();

      this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
      this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
      magnetic->resetMagneticField();
    }
  };

}}}}}

REGISTER_INTEGRATOR(
    Magnetic,Fixed,
    uammd::structured::Integrator::NVT::Magnetic::Fixed
)
