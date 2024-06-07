#ifndef __INTEGRATOR_BDHI_DOUBLY_PERIODIC_DPSOTKES__
#define __INTEGRATOR_BDHI_DOUBLY_PERIODIC_DPSOTKES__

#include <Integrator/BDHI/DoublyPeriodic/DPStokesSlab.cuh>

namespace uammd{
namespace structured{
namespace NVT{
namespace BDHIDoublyPeriodic{

	class DPStokes : public IntegratorBasicNVT{

		private:

			using BDHI = uammd::DPStokesSlab_ns::DPStokesIntegrator;

			std::unique_ptr<BDHI> bdhi;
			bool firstStep = true;

		public:

      DPStokes(std::shared_ptr<GlobalData>           gd,
               std::shared_ptr<ParticleGroup>        pg,
               DataEntry& data,
               std::string name):IntegratorBasicNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[DPStokes] This integrator can not handle more than one batch.");
				}

				BDHI::Parameters bdhiParameters;

				bdhiParameters.dt = this->dt;

				bdhiParameters.temperature = this->kBT;
				bdhiParameters.viscosity = data.getParameter<real>("viscosity");

				bdhiParameters.nx = data.getParameter<int>("nx");
				bdhiParameters.ny = data.getParameter<int>("ny");
				bdhiParameters.nz = data.getParameter<int>("nz",-1);

				bdhiParameters.Lx = data.getParameter<real>("Lx");
				bdhiParameters.Ly = data.getParameter<real>("Ly");
				bdhiParameters.H = data.getParameter<real>("H");

				bdhiParameters.w = data.getParameter<real>("w");
				bdhiParameters.beta  = data.getParameter<real>("beta",-1.0);
				bdhiParameters.alpha = data.getParameter<real>("alpha",-1.0);

				if(bdhiParameters.alpha > bdhiParameters.w*0.5){
					System::log<System::CRITICAL>("[DPStokes] alpha must be smaller than w/2");
				}

				bdhiParameters.w_d = data.getParameter<real>("w_d");
				bdhiParameters.beta_d  = data.getParameter<real>("beta_d",-1.0);
				bdhiParameters.alpha_d = data.getParameter<real>("alpha_d",-1.0);

				std::string mode = data.getParameter<std::string>("mode","none");
				if(mode == "none"){
					bdhiParameters.mode = uammd::DPStokesSlab_ns::WallMode::none;
				}else if(mode == "bottom"){
					bdhiParameters.mode = uammd::DPStokesSlab_ns::WallMode::bottom;
				}else if(mode == "slit"){
					bdhiParameters.mode = uammd::DPStokesSlab_ns::WallMode::slit;
				} else {
					System::log<System::CRITICAL>("[DPStokes] mode must be none, bottom or slit");
				}

				bdhiParameters.tolerance = data.getParameter<real>("tolerance",1e-7);

				System::log<System::MESSAGE>("[DPStokes] Temperature: %f",bdhiParameters.temperature);
				System::log<System::MESSAGE>("[DPStokes] Viscosity: %f",bdhiParameters.viscosity);

				System::log<System::MESSAGE>("[DPStokes] nx: %d",bdhiParameters.nx);
				System::log<System::MESSAGE>("[DPStokes] ny: %d",bdhiParameters.ny);
				System::log<System::MESSAGE>("[DPStokes] nz: %d",bdhiParameters.nz);

				System::log<System::MESSAGE>("[DPStokes] Lx: %f",bdhiParameters.Lx);
				System::log<System::MESSAGE>("[DPStokes] Ly: %f",bdhiParameters.Ly);
				System::log<System::MESSAGE>("[DPStokes] H: %f",bdhiParameters.H);

				System::log<System::MESSAGE>("[DPStokes] w: %f",bdhiParameters.w);
				System::log<System::MESSAGE>("[DPStokes] beta: %f",bdhiParameters.beta);
				System::log<System::MESSAGE>("[DPStokes] alpha: %f",bdhiParameters.alpha);

				System::log<System::MESSAGE>("[DPStokes] w_d: %f",bdhiParameters.w_d);
				System::log<System::MESSAGE>("[DPStokes] beta_d: %f",bdhiParameters.beta_d);
				System::log<System::MESSAGE>("[DPStokes] alpha_d: %f",bdhiParameters.alpha_d);

				System::log<System::MESSAGE>("[DPStokes] mode: %s",mode.c_str());

				System::log<System::MESSAGE>("[DPStokes] tolerance: %f",bdhiParameters.tolerance);

				bdhi = std::make_unique<BDHI>(pg->getParticleData(),
																			bdhiParameters);

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
