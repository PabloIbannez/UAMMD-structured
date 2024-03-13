//Template: Integrator

#ifndef __INTEGRATOR___CLASS_____SUBCLASS____
#define __INTEGRATOR___CLASS_____SUBCLASS____

namespace uammd{
namespace structured{
namespace __SUBTYPE__{
namespace __CLASS__{

	//For a standard integrator, you can use the following template.
	//Integrator basic only assumes that the integrator needs certain dt.
	//If you desired integrator needs define a dt (for a minimization, for example)
	//you can use the template below.

	//Time dependent integrator
	class __SUBCLASS__ : public IntegratorBasic{

		private:

			//Integrator basic is derived from uammd::Integrator.
			//Besides the uammd::Integrator attributes, it adds the following:
			//stream: the stream where the integrator will run. cudaStream_t
			//gd: a pointer to the global data. std::shared_ptr<GlobalData>
			//dt: time step. real

			//Declared here the variables that you need for the integrator
			//For example:

			real gamma;
			real alpha;

		public:

      __SUBCLASS__(std::shared_ptr<GlobalData>           gd,
                   std::shared_ptr<ParticleGroup>        pg,
                   DataEntry& data,
                   std::string name):IntegratorBasic(gd,pg,data,name){

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

	}

	//No time dependent integrator
	class __SUBCLASS__ : public Integrator{
			//This class derives from vanilla uammd::Integrator
		private:

			//To proper work with UAMMD-structured, probably you will need the following attributes:

			cudaStream_t stream;

			std::shared_ptr<GlobalData> gd;

			//Of course, you can add more attributes if you need them.
			//For example:

			real gamma;
			real alpha;

		public:

			__SUBCLASS__(std::shared_ptr<GlobalData>           gd,
									 std::shared_ptr<ParticleGroup>        pg,
									 DataEntry& data,
									 std::string name):Integrator(pg,name),gd(gd){

				//Set the stream

				stream = gd->getSystem()->getCudaStream();
				//stream is not used everywhere by default.
				//I recommend to use it when suming forces, for example.

				//Set up the integrator

				//You can read the parameters from data
				//For example:

				gamma = data.getParameter<real>("gamma");
				alpha = data.getParameter<real>("alpha",alphaDefault);

			}

			//You have to override the virtual function forwardTime.
			//Remember this function should also update the simulation step !!!
			//To do that, you can use the following line:
			//this->gd->setCurrentStep(this->gd->getCurrentStep()+1);

	}


}}}}

#endif
