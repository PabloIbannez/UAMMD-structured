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

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace DPD{

	class Verlet : public IntegratorBaseNVT{

		private:

			std::shared_ptr<uammd::VerletNVE> nve;
			bool firstStep = true;

		public:

      Verlet(std::shared_ptr<GlobalData>           gd,
             std::shared_ptr<ParticleGroup>        pg,
             DataEntry& data,
             std::string name):IntegratorBaseNVT(gd,pg,data,name){

				int batchNumber = GroupUtils::BatchGroupNumber(pg);
				if(batchNumber > 1){
					System::log<System::CRITICAL>("[Verlet] This integrator can not handle more than one batch.");
				}

				IntegratorUtils::generateVelocity(this->pg,
                                          this->kBT,
                                          this->gd->getSystem()->getSeed(),
                                          this->stream);

				typename uammd::VerletNVE::Parameters nveParams;

				nveParams.dt = this->dt;
				nveParams.initVelocities = false;

				nve = std::make_shared<uammd::VerletNVE>(this->pg, nveParams);

				//Set up DPD potential

				Potential::DPD::Parameters dpdParams;

				dpdParams.temperature = this->kBT;

				dpdParams.cutOff = data.getParameter<real>("cutOff");
				dpdParams.gamma  = data.getParameter<real>("gamma");

				dpdParams.A = data.getParameter<real>("alpha");

				dpdParams.dt = this->dt;

				System::log<System::MESSAGE>("[Verlet] cutOff: %f", dpdParams.cutOff);
				System::log<System::MESSAGE>("[Verlet] gamma: %f", dpdParams.gamma);
				System::log<System::MESSAGE>("[Verlet] alpha: %f", dpdParams.A);

				auto dpd = std::make_shared<Potential::DPD>(dpdParams);

				PairForces<Potential::DPD>::Parameters pairParams;

				pairParams.box = this->gd->getEnsemble()->getBox();

				auto dpdInteractor = std::make_shared<PairForces<Potential::DPD>>(this->pg->getParticleData(),
																																					pairParams,
																																					dpd);

				nve->addInteractor(dpdInteractor);

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
    DPD,Verlet,
    uammd::structured::Integrator::NVT::DPD::Verlet
)
