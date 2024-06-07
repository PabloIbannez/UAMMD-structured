#ifndef __INTEGRATOR_BDHI_TRIPLY_PERIODIC_FIB__
#define __INTEGRATOR_BDHI_TRIPLY_PERIODIC_FIB__

#include"Integrator/BDHI/FIB.cuh"

namespace uammd{
namespace structured{
namespace NVT{
namespace BDHITriplyPeriodic{

	class FluctuatingImmersedBoundary : public IntegratorBasicNVT{

		private:

			using BDHI = BDHI::FIB;

			std::unique_ptr<BDHI> bdhi;
			bool firstStep = true;

		public:

      FluctuatingImmersedBoundary(std::shared_ptr<GlobalData>           gd,
               					          std::shared_ptr<ParticleGroup>        pg,
               					          DataEntry& data,
               					          std::string name):IntegratorBasicNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[BDHI] This integrator can not handle more than one batch.");
				}

				BDHI::Parameters bdhiParameters;

				bdhiParameters.dt = this->dt;

				bdhiParameters.temperature = this->kBT;
				bdhiParameters.viscosity = data.getParameter<real>("viscosity");

				bdhiParameters.hydrodynamicRadius = data.getParameter<real>("hydrodynamicRadius",1.0);

				bdhiParameters.tolerance = data.getParameter<real>("tolerance",1e-3);
				bdhiParameters.scheme    = uammd::BDHI::FIB::Scheme::IMPROVED_MIDPOINT;

				bdhiParameters.box = gd->getEnsemble()->getBox();

				System::log<System::MESSAGE>("[BDHI] Viscosity: %f",bdhiParameters.viscosity);
				System::log<System::MESSAGE>("[BDHI] Hydrodynamic radius: %f",bdhiParameters.hydrodynamicRadius);
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

}}}}

#endif
