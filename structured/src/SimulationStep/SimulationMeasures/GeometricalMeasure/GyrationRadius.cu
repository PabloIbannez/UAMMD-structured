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

class GyrationRadius: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

    public:

        GyrationRadius(std::shared_ptr<ParticleGroup>             pg,
                       std::shared_ptr<IntegratorManager> integrator,
                       std::shared_ptr<ForceField>                ff,
                       DataEntry& data,
                       std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            //Read parameters

            outputFilePath = data.getParameter<std::string>("outputFilePath");
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step GyrationRadius" << std::endl;
            }
        }

        void applyStep(ullint step, cudaStream_t st) override{

            Box box = gd->getEnsemble()->getBox();

            outputFile << step << " ";

            real3 centroid = Measures::centroidPos(pg);
            real  Rg       = Measures::gyrationRadius(pg,centroid,box);

            outputFile << Rg << std::endl;
        }
};

}}}}

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,GyrationRadius,
    uammd::structured::SimulationStep::SimulationMeasures::GyrationRadius
)
