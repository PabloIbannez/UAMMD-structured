#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace Special{
namespace None{

	class  None: public IntegratorBase{

		public:

      None(std::shared_ptr<GlobalData>    gd,
           std::shared_ptr<ParticleGroup> pg,
           DataEntry& data,
           std::string name):IntegratorBase(gd,pg,data,name){
			}

			void forwardTime() override {
        this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
			}

	};

}}}}}

REGISTER_INTEGRATOR(
    None,None,
    uammd::structured::Integrator::Special::None::None
)
