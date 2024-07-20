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

class Height: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        int N;
        int particleNumberAverage;

    public:

        Height(std::shared_ptr<ParticleGroup>             pg,
               std::shared_ptr<IntegratorManager> integrator,
               std::shared_ptr<ForceField>                ff,
               DataEntry& data,
               std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            outputFilePath = data.getParameter<std::string>("outputFilePath");

            particleNumberAverage = data.getParameter<int>("particleNumberAverage",1);
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step height ..." << std::endl;
            }

            N = pg->getNumberParticles();
        }

        void applyStep(ullint step, cudaStream_t st) override{

            outputFile << step << " ";

            auto pos = pd->getPos(access::location::cpu, access::mode::read);

            auto posIterator = pg->getPropertyIterator(pd->getPos(access::location::cpu, access::mode::read).begin(),
                                                       access::location::cpu);

            //Sort positions by z coordinate, from highest to lowest. Not in place.
            std::vector<real4> sortedPos(posIterator, posIterator + N);
            std::sort(sortedPos.begin(), sortedPos.end(), [](real4 a, real4 b){return a.z > b.z;});

            //Compute average height of the particleNumberAverage particles with highest z coordinate.
            real height = 0;
            for(int i = 0; i < particleNumberAverage; i++){
                height += sortedPos[i].z;
            }

            //Write average height to file.
            outputFile << height/particleNumberAverage << std::endl;

        }
};

}}}}

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,Height,
    uammd::structured::SimulationStep::SimulationMeasures::Height
)
