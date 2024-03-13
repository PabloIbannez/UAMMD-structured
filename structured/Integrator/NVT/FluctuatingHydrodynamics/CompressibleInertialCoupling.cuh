#ifndef __INTEGRATOR_FLUCTUATING_HYDRODYNAMICS__COMPRESSIBLE_INERTIAL_COUPLING__
#define __INTEGRATOR_FLUCTUATING_HYDRODYNAMICS__COMPRESSIBLE_INERTIAL_COUPLING__

#include"Integrator/Hydro/ICM_Compressible.cuh"

namespace uammd{
namespace structured{
namespace NVT{
namespace FluctuatingHydrodynamics{

	class CompressibleInertialCoupling : public IntegratorBasicNVT{

		private:

			using ICM = Hydro::ICM_Compressible;

			std::unique_ptr<ICM> icm;
			bool firstStep = true;

		public:

      CompressibleInertialCoupling(std::shared_ptr<GlobalData>           gd,
               					  				 std::shared_ptr<ParticleGroup>        pg,
               					  				 DataEntry& data,
               					  				 std::string name):IntegratorBasicNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[ICM] This integrator can not handle more than one batch.");
				}

				// Check particle group is composed of all particles
				if(pg->getNumberParticles() != pg->getParticleData()->getNumParticles()){
						System::log<System::CRITICAL>("[ICM] Only groups composed of all particles are allowed.");
				}

				ICM::Parameters icmParameters;

				icmParameters.dt = this->dt;
				icmParameters.temperature = this->kBT;

				icmParameters.shearViscosity = data.getParameter<real>("shearViscosity");
				icmParameters.bulkViscosity  = data.getParameter<real>("bulkViscosity");

				icmParameters.speedOfSound = data.getParameter<real>("speedOfSound");
				real density               = data.getParameter<real>("density");

				icmParameters.hydrodynamicRadius = data.getParameter<real>("hydrodynamicRadius");

				icmParameters.boxSize = gd->getEnsemble()->getBox().boxSize;

				int3 cellDim;

				cellDim.x = data.getParameter<int>("cellDimX");
				cellDim.y = data.getParameter<int>("cellDimY");
				cellDim.z = data.getParameter<int>("cellDimZ");

				icmParameters.cellDim = cellDim;

				icmParameters.seed = gd->getSystem()->getSeed();

				std::string init = data.getParameter<std::string>("initialCondition","none");

				//You can add more initial conditions here
				if(init == "none"){
					icmParameters.initialDensity   = [=](real3 position){return density;};
					icmParameters.initialVelocityX = [](real3 position){return 0.0;};
					icmParameters.initialVelocityY = [](real3 position){return 0.0;};
					icmParameters.initialVelocityZ = [](real3 position){return 0.0;};
				} else {
					System::log<System::CRITICAL>("[ICM] Initial condition not recognized. Available options are: none");
				}

				icm = std::make_unique<ICM>(pg->getParticleData(),
																		icmParameters);

			}

			void forwardTime() override {

				if(firstStep){
					//Load all interactors into icm
					for(auto& interactor : this->getInteractors()){
						icm->addInteractor(interactor);
					}

					//Load all updatables into icm
					for(auto& updatable : this->getUpdatables()){
						icm->addUpdatable(updatable);
					}

					firstStep = false;
				}

				icm->forwardTime();
        this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
			}

	};

}}}}

#endif
