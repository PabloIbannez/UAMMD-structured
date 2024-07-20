#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

#include "Integrator/VerletNVE.cuh"
#include "Interactor/SPH.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace SPH{

	class Verlet : public IntegratorBaseNVT{

		private:

			std::shared_ptr<uammd::VerletNVE> nve;
			bool firstStep = true;

		public:

      Verlet(std::shared_ptr<GlobalData>           gd,
             std::shared_ptr<ParticleGroup>        pg,
             DataEntry& data,
             std::string name):IntegratorBaseNVT(gd,pg,data,name){

				IntegratorUtils::generateVelocity(this->pg,
                                          this->kBT,
                                          this->gd->getSystem()->getSeed(),
                                          this->stream);

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[Verlet] This integrator can not handle more than one batch.");
				}

				typename uammd::VerletNVE::Parameters nveParams;

				nveParams.dt = this->dt;
				nveParams.initVelocities = false;

				nve = std::make_shared<uammd::VerletNVE>(this->pg, nveParams);

				//Set up SPH potential

				uammd::SPH::Parameters sphParams;

				sphParams.box = this->gd->getEnsemble()->getBox();

				sphParams.support = data.getParameter<real>("support");

				sphParams.viscosity = data.getParameter<real>("viscosity");

				sphParams.gasStiffness = data.getParameter<real>("gasStiffness");
				sphParams.restDensity  = data.getParameter<real>("restDensity");

				System::log<System::MESSAGE>("[Verlet] support: %f", sphParams.support);
				System::log<System::MESSAGE>("[Verlet] viscosity: %f", sphParams.viscosity);
				System::log<System::MESSAGE>("[Verlet] gasStiffness: %f", sphParams.gasStiffness);
				System::log<System::MESSAGE>("[Verlet] restDensity: %f", sphParams.restDensity);

				auto sph = std::make_shared<uammd::SPH>(this->pg,sphParams);

				nve->addInteractor(sph);

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
    SPH,Verlet,
    uammd::structured::Integrator::NVT::SPH::Verlet
)
