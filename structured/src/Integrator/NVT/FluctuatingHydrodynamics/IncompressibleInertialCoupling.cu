#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

#include"Integrator/VerletNVE.cuh"
#include"Interactor/PairForces.cuh"
#include"Interactor/Potential/DPD.cuh"

#include"Integrator/Hydro/ICM.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace FluctuatingHydrodynamics{

	class IncompressibleInertialCoupling : public IntegratorBaseNVT{

		private:

			using ICM = Hydro::ICM;

			std::unique_ptr<ICM> icm;
			bool firstStep = true;

		public:

      IncompressibleInertialCoupling(std::shared_ptr<GlobalData>           gd,
                 					  				 std::shared_ptr<ParticleGroup>        pg,
                 					  				 DataEntry& data,
                 					  				 std::string name):IntegratorBaseNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[ICM] This integrator can not handle more than one batch.");
				}

				ICM::Parameters icmParameters;

				icmParameters.dt = this->dt;
				icmParameters.temperature = this->kBT;

				icmParameters.viscosity = data.getParameter<real>("viscosity");
				icmParameters.density   = data.getParameter<real>("density");

				icmParameters.hydrodynamicRadius = data.getParameter<real>("hydrodynamicRadius");

				icmParameters.box = gd->getEnsemble()->getBox();

				if(data.isParameterAdded("cellDimX") or
					 data.isParameterAdded("cellDimY") or
					 data.isParameterAdded("cellDimZ")){

					int3 cellDim;

					cellDim.x = data.getParameter<int>("cellDimX");
					cellDim.y = data.getParameter<int>("cellDimY");
					cellDim.z = data.getParameter<int>("cellDimZ");

					icmParameters.cells = cellDim;
				}

				icmParameters.sumThermalDrift     = data.getParameter<bool>("sumThermalDrift",false);
				icmParameters.removeTotalMomentum = data.getParameter<bool>("removeTotalMomentum",true);

				icm = std::make_unique<ICM>(pg,
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

}}}}}

REGISTER_INTEGRATOR(
    FluctuatingHydrodynamics,IncompressibleInertialCoupling,
    uammd::structured::Integrator::NVT::FluctuatingHydrodynamics::IncompressibleInertialCoupling
)
