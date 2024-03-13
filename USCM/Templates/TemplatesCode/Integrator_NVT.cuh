//Template: Integrator NVT

#ifndef __INTEGRATOR___CLASS_____SUBCLASS____
#define __INTEGRATOR___CLASS_____SUBCLASS____

namespace uammd{
namespace structured{
namespace __SUBTYPE__{
namespace __CLASS__{

	//IntegratorBasicNVT is simply an IntegratorBasic with two additional parameters
	//that are the temperature and kBT. Since temperature is also stored by GlobalData
	//it set ups the integrator taking this into account. If temperature is given in
	//parameters, it is used instead of the one in GlobalData otherwise the one in
	//GlobalData is used.
	//kBT is computed after temperature is set up. kBT depends of the units used.
	//kBT = gd->getUnits()->getConstant("KBOLTZ")*temperature

	class __SUBCLASS__ : public IntegratorBasicNVT{

		private:

			//You can declare additional parameters here.
			//For example:
			real alpha;
			real gamma;

		public:

      __SUBCLASS__(std::shared_ptr<GlobalData>           gd,
                   std::shared_ptr<ParticleGroup>        pg,
                   DataEntry& data,
                   std::string name):IntegratorBasicNVT(gd,pg,data,name){

				//Set up the integrator

				//You can read the parameters from data
				//For example:

				gamma = data.getParameter<real>("gamma");
				alpha = data.getParameter<real>("alpha",alphaDefault);

			}

			//IntegratorBasic adds the following methods (plus the ones from uammd::Integrator),
			//feel free to use them or not:
			//void resetForce(): sets the force to zero
			//void resetEnergy(): sets the energy to zero
			//void updateForce(): updates the force, suming all the interactors
			//void updateEnergy(): updates the energy, suming all the interactors
			//ALL THE PREVIOUS METHODS USE THE STREAM stream !!!

			//Now you have to override the virtual function forwardTime

			void forwardTime() override {

				//Perform the integration step

				//...
				//...
				//...

				//Remember integrators has to update the simulation step !!!
				//This can be done with the following line
        this->gd->setCurrentStep(this->gd->getCurrentStep()+1);
			}

	};

}}}}

#endif
