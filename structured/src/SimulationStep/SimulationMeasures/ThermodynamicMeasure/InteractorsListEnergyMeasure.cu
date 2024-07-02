#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class InteractorsListEnergyMeasure: public SimulationStepBase_EnergyForceTorque{

        std::string    outputFilePath;
        std::ofstream  outputFile;

        std::vector<std::string>                                          selectedInteractorNamesList;
        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> selectedInteractors;

    public:

        InteractorsListEnergyMeasure(std::shared_ptr<ParticleGroup>             pg,
                                     std::shared_ptr<IntegratorManager> integrator,
                                     std::shared_ptr<ForceField>                ff,
                                     DataEntry& data,
                                     std::string name):SimulationStepBase_EnergyForceTorque(pg,
                                                                                            integrator,ff,
                                                                                            data,name){

            //Read parameters
            outputFilePath              = data.getParameter<std::string>("outputFilePath");
            if(data.isParameterAdded("interactorsList")){
                selectedInteractorNamesList = data.getParameter<std::vector<std::string>>("interactorsList");
            } else {
                selectedInteractorNamesList = std::vector<std::string>();
            }
        }


        void init(cudaStream_t st) override {

            //

            std::map<std::string,std::shared_ptr<typename uammd::Interactor>> allInteractors = this->topology->getInteractors();

            std::vector<std::string> allInteractorsNames;
            for(auto& interactor:allInteractors){
                allInteractorsNames.push_back(interactor.first);
            }

            //If no interactors are selected, use all of them
            if(selectedInteractorNamesList.size()==0){
                selectedInteractorNamesList = allInteractorsNames;
            }

            //Check all selected interactors are in the list of all interactors
            for(std::string& interactorName : selectedInteractorNamesList){
                if(std::find(allInteractorsNames.begin(),allInteractorsNames.end(),interactorName) == allInteractorsNames.end()){
                    System::log<System::CRITICAL>("[InteractorsListEnergyMeasure] Selected interactor %s is not in present in force field",
                                                  interactorName.c_str());
                }
            }

            //Get the interactors
            for(std::string& interactorName : selectedInteractorNamesList){
                selectedInteractors[interactorName] = allInteractors[interactorName];
            }

            // Write header

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << std::setw(12) << "# Step";

                for(auto& interactor: selectedInteractors){
                    outputFile << std::setw(12) << interactor.first << " ";
                }
                outputFile << std::endl;
            }
        }

        void applyStep(ullint step, cudaStream_t st) override{

            //Reset energy, force
            this->setZero(st);

            //Write current step to output file
            outputFile << std::setw(12) << step << " ";

            //Compute energy
				    uammd::Interactor::Computables comp;
				    comp.energy = true;

            //Sum energy
    				for(auto& interactor: selectedInteractors){
                this->setZero(st);
                interactor.second->sum(comp,st);

                real totalEnergy = Measures::totalPotentialEnergy(pg);

                outputFile << std::setw(12) << totalEnergy << " ";
            }
            outputFile << std::endl;
        }
};

}}}}

REGISTER_SIMULATION_STEP(
    ThermodynamicMeasure,InteractorsListEnergyMeasure,
    uammd::structured::SimulationStep::SimulationMeasures::InteractorsListEnergyMeasure
)
