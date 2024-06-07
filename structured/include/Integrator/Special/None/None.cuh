#ifndef __INTEGRATOR_NONE__
#define __INTEGRATOR_NONE__

namespace uammd{
namespace structured{
namespace Special{
namespace None{

	class  None: public IntegratorBasic{

		public:

      None(std::shared_ptr<GlobalData>    gd,
           std::shared_ptr<ParticleGroup> pg,
           DataEntry& data,
           std::string name):IntegratorBasic(gd,pg,data,name){
			}

			void forwardTime() override {
        this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
			}

	};

}}}}

#endif
