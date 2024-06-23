#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

#include "Integrator/VerletNVE.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVE{
namespace Verlet{

	class VelocityVerlet : public IntegratorBaseNVT{

		private:

			std::unique_ptr<uammd::VerletNVE> nve;
			bool firstStep = true;

		public:

      VelocityVerlet(std::shared_ptr<GlobalData>           gd,
                     std::shared_ptr<ParticleGroup>        pg,
                     DataEntry& data,
                     std::string name):IntegratorBaseNVT(gd,pg,data,name){

				IntegratorUtils::generateVelocity(this->pg,
                                                  this->kBT,
                                                  this->gd->getSystem()->getSeed(),
                                                  this->stream);

				typename uammd::VerletNVE::Parameters nveParams;

				nveParams.dt = this->dt;
				nveParams.initVelocities = false;

				nve = std::make_unique<uammd::VerletNVE>(this->pg, nveParams);

			}

			void forwardTime() override {

				if(firstStep){
					//Load all interactors into nve
					for(auto& interactor : this->getInteractors()){
						nve->addInteractor(interactor);
					}

					//Load all updatables into nve
					for(auto& updatable : this->getUpdatables()){
						nve->addUpdatable(updatable);
					}

					firstStep = false;
				}

				nve->forwardTime();
        this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
			}

	};

}}}}}

REGISTER_INTEGRATOR(
    Verlet,VelocityVerlet,
    uammd::structured::Integrator::NVE::Verlet::VelocityVerlet
)
