#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

#include "Integrator/BDHI/BDHI_EulerMaruyama.cuh"
#include "Integrator/BDHI/BDHI_Lanczos.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace BDHIOpenBoundary{

	class Lanczos : public IntegratorBaseNVT{

		private:

			using BDHI = BDHI::EulerMaruyama<BDHI::Lanczos>;

			std::unique_ptr<BDHI> bdhi;
			bool firstStep = true;

		public:

      Lanczos(std::shared_ptr<GlobalData>           gd,
              std::shared_ptr<ParticleGroup>        pg,
              DataEntry& data,
              std::string name):IntegratorBaseNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[BDHI] This integrator can not handle more than one batch.");
				}

				BDHI::Parameters bdhiParameters;

				bdhiParameters.dt = this->dt;

				bdhiParameters.temperature = this->kBT;
				bdhiParameters.viscosity = data.getParameter<real>("viscosity");

				bdhiParameters.hydrodynamicRadius = data.getParameter<real>("hydrodynamicRadius",-1.0);

				bdhiParameters.tolerance = data.getParameter<real>("tolerance",1e-3);

				System::log<System::MESSAGE>("[BDHI] Viscosity: ",bdhiParameters.viscosity);

				if(bdhiParameters.hydrodynamicRadius < 0.0){
					System::log<System::MESSAGE>("[BDHI] Hydrodynamic radius not set, using particle radius");
				} else {
					System::log<System::MESSAGE>("[BDHI] Hydrodynamic radius: %f",bdhiParameters.hydrodynamicRadius);
				}

				System::log<System::MESSAGE>("[BDHI] Tolerance: %f",bdhiParameters.tolerance);

				bdhi = std::make_unique<BDHI>(pg,bdhiParameters);

			}

			void forwardTime() override {

				if(firstStep){
					//Load all interactors into bdhi
					for(auto& interactor : this->getInteractors()){
						bdhi->addInteractor(interactor);
					}

					//Load all updatables into bdhi
					for(auto& updatable : this->getUpdatables()){
						bdhi->addUpdatable(updatable);
					}

					firstStep = false;
				}

				bdhi->forwardTime();
        this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
			}

	};

}}}}}

REGISTER_INTEGRATOR(
    BDHIOpenBoundary,Lanczos,
    uammd::structured::Integrator::NVT::BDHIOpenBoundary::Lanczos
)
