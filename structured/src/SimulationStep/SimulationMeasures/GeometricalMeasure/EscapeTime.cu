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

class insideSystemChecker{
        private:
                std::vector<real3> normalVector;
                std::vector<real3> independentVector;
        public:
                insideSystemChecker(std::vector<real3> normalVector, std::vector<real3> independentVector):normalVector(normalVector),independentVector(independentVector){}
                bool isInsideSystem(real3 pos){ //Check if the particle is below all the planes contained in normalVector and independentVector. By computing the dot product of normalVector * (pos-independentVector)
                    bool isInside = true;
                    for (int plane = 0; plane < normalVector.size(); plane++){
                        real3 particlePositionMinusIndependentVector = pos - independentVector[plane];
                        real dotProduct = dot(normalVector[plane],particlePositionMinusIndependentVector);

                        if(dotProduct < 0){
                            isInside = false;
                        }
                    }
                    return isInside;
                }
};

class EscapeTime: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        std::shared_ptr<insideSystemChecker> checker;
        std::map<int,bool> isParticleInsideSystem;

    public:

        EscapeTime(std::shared_ptr<ParticleGroup>         pg,
               std::shared_ptr<IntegratorManager> integrator,
               std::shared_ptr<ForceField>                ff,
               DataEntry& data,
               std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            //Read parameters from input file.
            outputFilePath = data.getParameter<std::string>("outputFilePath");

            //Read Data from input file.
            auto normalVectorData      = data.getData<real3>("normalVector");
            auto independentVectorData = data.getData<real3>("independentVector");

            std::vector<real3> normalVector;
            std::vector<real3> independentVector;

            for (int i = 0; i < normalVectorData.size(); i++){
                normalVector.push_back({normalVectorData[i].x,normalVectorData[i].y,normalVectorData[i].z});
                independentVector.push_back({independentVectorData[i].x,independentVectorData[i].y,independentVectorData[i].z});
            }

            checker = std::make_shared<insideSystemChecker>(normalVector,independentVector);

            auto pos = pd->getPos(access::location::cpu, access::mode::read);
            const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
            for(int i = 0; i < pos.size(); i++){
                int index = sortedIndex[i];
                real3 particlePosition = make_real3(pos[index]);
                isParticleInsideSystem[index] = checker->isInsideSystem(particlePosition);
                if (!isParticleInsideSystem[index]){
                    outputFile << 0 << " " << index << std::endl;
                }
            }
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);
        }

        void applyStep(ullint step, cudaStream_t st) override{

            auto pos = pd->getPos(access::location::cpu, access::mode::read);
            const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
            for(int i = 0; i < pos.size(); i++){
                int index = sortedIndex[i];
                real3 particlePosition = make_real3(pos[index]);
                bool isInsideSystem = checker->isInsideSystem(particlePosition);

                if(isParticleInsideSystem[index] && !isInsideSystem){
                    outputFile << step << " " << index << std::endl;
                    isParticleInsideSystem[index] = false;
                }
            }

        }
};

}}}}

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,EscapeTime,
    uammd::structured::SimulationStep::SimulationMeasures::EscapeTime
)
