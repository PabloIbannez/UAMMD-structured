//Template: Integrator

#ifndef __INTEGRATOR___CLASS_____SUBCLASS____
#define __INTEGRATOR___CLASS_____SUBCLASS____

namespace uammd{
namespace structured{
namespace __TYPE__{
namespace __SUBTYPE__{

	//Simulation step template

	class __SUBCLASS__ : public SimulationStepBase{

		private:

			//It is common for a simulation step to write to a file.
			//We can declare the file here, and open it in the constructor.

			std::string   outputFilePath;
			std::ofstream outputFile;

			//If additional data is needed, we can declare it here.
			//For example:

			real data1;
			real data2;

		public:

			//Simulation step has all the information about the simulation.
			//ParticleGroup, pg: Contains the particles of the current group.
			//Take into accunt that we can ask pg for ParticleData, with pg->getParticleData().
			//IntegratorManager, integrator: Contains all the integrators of the simulation.
			//ForceField::ForceFieldBase, forcefield: Contains all the interactor of the simulation.

			//DataEntry, is a common DataEntry type. The following parameters are handled by the SimulationStepBase:
			// - startStep: The first step where the simulation step will be executed.
			// - endStep: The last step where the simulation step will be executed.
			// - intervalStep: The interval between executions.

			//If you want to access them use getters, for example:
			// - this->getStartStep()
			// - this->getEndStep()
			// - this->getIntervalStep()

			//You have also the following methods:
			// - this->getName(): Returns the name of the simulation step.
			// - this->getLastStepApplied(): Returns the last step where the simulation step was executed.

			//Additonal classes are also available (they are requested to pg, integrator and forcefield):
			//topology: Contains the topology of the simulation.
			//sys: Contains the system of the simulation.
			//gd: Contains the global data of the simulation.
			//pd: Contains the particle data of the simulation.

			//You can use this->topology, this->sys, this->gd and this->pd to access them.

      __SUBCLASS__(std::shared_ptr<ParticleGroup>  pg,
                   std::shared_ptr<IntegratorManager> integrator,
                   std::shared_ptr<ForceField::ForceFieldBase>    ff,
                   DataEntry& data,
                   std::string name):SimulationStepBase(pg,integrator,ff,data,name){

				//Set the output file path.
				//We can use the DataEntry to get the output file path.
				outputFilePath = data.getParameter<std::string>("outputFilePath");
				//Is better not to open the file here, but in init function.
				//At init function we can ensure pg, integrator and forcefield are initialized.

				//We can also get additional data from the DataEntry.
				//For example:
				data1 = data.getParameter<real>("data1");

			}

			//We have to override the init function.
			void init(cudaStream_t st) override{

				//It is recommended to open the file here.
				//It is also recommended to check if the file exists.
				//It can be done with the function Backup::openFile.

				bool isFileOpen = Backup::openFile(this->sys, outputFilePath, outputFile);

				//This function opens the file.
				//If the file exists, it will be opened in append mode and isFileOpen will be true.
				//If the file does not exist, it will be created and isFileOpen will be false.
				//If there is an error, the function will throw an exception.

				//If the file did not exist, we can write the header here.
				if(!isFileOpen){
					outputFile << "Step Measured1 Measured2 ..." << std::endl;
				}

				//We can also init here additional variables.
				//Now we can be sure that pg, integrator and forcefield are initialized.

				//For example:
				data2 = Measures::totalMass(pg);

			}

			//We have to override the applyStep function.

			void applyStep(ullint step, cudaStream_t st) override{

				//Now is time to compute the data we want to write to the file.
				//For example:
				real com = Measures::centerOfMassPos(pg,st);
				//Take into account thant if you want use the gpu, you should use the stream st.

				//Operate with the data.
				com = com + data1 + data2;

				//Write the data to the file.

				outputFile << step << " " << com << std::endl;

			}

	};

}}}}

#endif
