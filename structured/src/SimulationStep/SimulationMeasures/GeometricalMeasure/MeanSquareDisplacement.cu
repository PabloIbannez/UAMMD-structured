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

class MeanSquareDisplacement: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        std::map<int,real3> referencePositions;
        bool firstStep = true;
        //Only in the first step, the reference positions are stored

    public:

        MeanSquareDisplacement(std::shared_ptr<ParticleGroup>  pg,
                 std::shared_ptr<IntegratorManager> integrator,
                 std::shared_ptr<ForceField>                ff,
                 DataEntry& data,
                 std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            outputFilePath = data.getParameter<std::string>("outputFilePath");
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step MSD" << std::endl;
            }

        }

        void applyStep(ullint step, cudaStream_t st) override{

            std::map<int,int> id_index;

            {
                auto id   = pd->getId(access::location::cpu,
                                      access::mode::read);

                auto groupIndex  = pg->getIndexIterator(access::location::cpu);
                auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                fori(0,pg->getNumberParticles()){
                    int id_   = id[groupIndex[i]];
                    int index = sortedIndex[id_];

                    id_index[id_]=index;
                }
            }

            if(firstStep){

                auto pos = pd->getPos(access::location::cpu,
                                      access::mode::read);

                for(auto& id_i : id_index){
                    int id     = id_i.first;
                    int index  = id_i.second;

                    referencePositions[id] = make_real3(pos[index]);
                }

                firstStep = false;
            }

            //Write current step to output file
            outputFile << step << " ";

            //Compute the mean square displacement
            real msd = 0;

            auto id  = pd->getId(access::location::cpu,
                                  access::mode::read);

            auto currentPositions = this->pd->getPos(access::location::cpu,access::mode::read);

            auto groupIndex       = pg->getIndexIterator(access::location::cpu);

            fori(0,pg->getNumberParticles()){
                int id_   = id[groupIndex[i]];
                int index = id_index[id_];

                real3 dr = make_real3(currentPositions[index])  - referencePositions[id_];

                msd += dot(dr,dr);
            }

            msd /= pg->getNumberParticles();

            outputFile << msd << std::endl;
        }
};

}}}}

REGISTER_SIMULATION_STEP(
    GeometricalMeasure,MeanSquareDisplacement,
    uammd::structured::SimulationStep::SimulationMeasures::MeanSquareDisplacement
)
