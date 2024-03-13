#ifndef __SIMULATION_STEP_PARTICLE_LIST_DISTANCES_MEASURES__
#define __SIMULATION_STEP_PARTICLE_LIST_DISTANCES_MEASURES__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

    class DistancesMeasure: public SimulationStepBase{

            std::string   outputFilePath;
            std::ofstream outputFile;

            std::vector<int2> ids2;

        public:

            DistancesMeasure(std::shared_ptr<ParticleGroup>             pg,
                             std::shared_ptr<IntegratorManager> integrator,
                             std::shared_ptr<ForceField>                ff,
                             DataEntry& data,
                             std::string name):SimulationStepBase(pg,integrator,ff,data,name){

                //Read parameters

                outputFilePath = data.getParameter<std::string>("outputFilePath");

                //////////////////////////////////////////////

                auto id_i = data.getData<int>("id_i");
                auto id_j = data.getData<int>("id_j");

                for(int index=0;index<data.getDataSize();index++){
                    ids2.push_back({id_i[index],id_j[index]});
                }
            }

            void init(cudaStream_t st) override{
                bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);
            }

            void applyStep(ullint step, cudaStream_t st) override{

                Box box = gd->getEnsemble()->getBox();

                auto pos = pd->getPos(access::location::cpu,
                                      access::mode::read);

                const int *sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                for(auto ids : ids2){
                    int i = sortedIndex[ids.x];
                    int j = sortedIndex[ids.y];

                    real3 posi = make_real3(pos[i]);
                    real3 posj = make_real3(pos[j]);

                    real r = length(box.apply_pbc(posj-posi));

                    outputFile << r << " ";
                }

                outputFile << std::endl;
            }
    };

}}}}

#endif
