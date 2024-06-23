#ifndef __INTEGRATOR_PARTICLESLIST_MEASURE_POTENTIAL_MEASURE__
#define __INTEGRATOR_PARTICLESLIST_MEASURE_POTENTIAL_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

	class PotentialMeasure : public SimulationStepBase_EnergyForceTorque{

		private:

			std::string   outputFilePath;
			std::ofstream outputFile;

			std::vector<int> ids;

			std::vector<std::string> interactorNames;

			std::map<int,std::map<std::string,real>>  energyData;
			std::map<int,std::map<std::string,real4>> forceData;
			std::map<int,std::map<std::string,real4>> torqueData;

		public:

      PotentialMeasure(std::shared_ptr<ParticleGroup>             pg,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField>    					  ff,
                       DataEntry& data,
                       std::string name):SimulationStepBase_EnergyForceTorque(pg,
												 																											integrator,ff,
																																							data,name){

				outputFilePath = data.getParameter<std::string>("outputFilePath");

				//Read ids of selected particles
				ids = data.getData<int>("id");

			}

			void init(cudaStream_t st) override{

				//Get interactor names. this->topology->getInteractors() returns a map, key: interactor name, value: interactor shared pointer
				for(auto& interactor : this->topology->getInteractors()){
					interactorNames.push_back(interactor.first);
				}

				bool isFileEmpty = Backup::openFile(this->sys, outputFilePath, outputFile);

				if(isFileEmpty){
					//Write header
					outputFile << std::setw(24) << "# id ";
					for(std::string& interactorName : interactorNames){
						outputFile << std::setw(24) << interactorName + "Energy ";

						outputFile << std::setw(24) << interactorName + "ForceX ";
						outputFile << std::setw(24) << interactorName + "ForceY ";
						outputFile << std::setw(24) << interactorName + "ForceZ ";
						outputFile << std::setw(24) << interactorName + "ForceMod ";

						outputFile << std::setw(24) << interactorName + "TorqueX ";
						outputFile << std::setw(24) << interactorName + "TorqueY ";
						outputFile << std::setw(24) << interactorName + "TorqueZ ";
						outputFile << std::setw(24) << interactorName + "TorqueMod ";

					}
					outputFile << std::endl;
				}

				//Set the precision of the output file
				if constexpr (std::is_same<real, float>::value) {
					outputFile << std::fixed << std::setprecision(std::numeric_limits<float>::digits10);
				} else if constexpr (std::is_same<real, double>::value) {
					outputFile << std::fixed << std::setprecision(std::numeric_limits<double>::digits10);
				}

			}

			void applyStep(ullint step, cudaStream_t st) override {

				outputFile << "# " << step << "\n";

				auto id2index = pd->getIdOrderedIndices(access::location::cpu);

				uammd::Interactor::Computables comp;
				comp.energy = true;
				comp.force  = true;

				for(std::string& interactorName : interactorNames){
					this->setZero(st);
					auto interactor = this->topology->getInteractor(interactorName);
					interactor->sum(comp, st);

					for(int i = 0; i < ids.size(); i++){
						int index = id2index[ids[i]];

						auto energy = this->pd->getEnergy(access::location::cpu, access::mode::read);
						energyData[i][interactorName] = energy[index];

						auto force = this->pd->getForce(access::location::cpu, access::mode::read);
						forceData[i][interactorName] = force[index];

						auto torque = this->pd->getTorque(access::location::cpu, access::mode::read);
						torqueData[i][interactorName] = torque[index];
					}
				}

				//Write data

				for(int i = 0; i < ids.size(); i++){
					outputFile << std::setw(24) << ids[i];

					for(std::string& interactorName : interactorNames){
						outputFile << std::setw(24) << energyData[i][interactorName];

						real3 force = make_real3(forceData[i][interactorName]);
						outputFile << std::setw(24) << force.x;
						outputFile << std::setw(24) << force.y;
						outputFile << std::setw(24) << force.z;
						outputFile << std::setw(24) << length(force);

						real3 torque = make_real3(torqueData[i][interactorName]);
						outputFile << std::setw(24) << torque.x;
						outputFile << std::setw(24) << torque.y;
						outputFile << std::setw(24) << torque.z;
						outputFile << std::setw(24) << length(torque);
					}

					outputFile << std::endl;
				}

			}

	};

}}}}

#endif
