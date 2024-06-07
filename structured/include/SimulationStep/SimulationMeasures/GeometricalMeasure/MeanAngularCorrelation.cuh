#ifndef __SIMULATION_MEASURES_MEAN_ANGULAR_CORRELATION__
#define __SIMULATION_MEASURES_MEAN_ANGULAR_CORRELATION__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class MeanAngularCorrelation: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        std::map<int,real4> referenceDirections;
        bool firstStep = true;
        //Only in the first step, the reference directions are stored

    public:

        MeanAngularCorrelation(std::shared_ptr<ParticleGroup>             pg,
                               std::shared_ptr<IntegratorManager> integrator,
                               std::shared_ptr<ForceField>                ff,
                               DataEntry& data,
                               std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            outputFilePath = data.getParameter<std::string>("outputFilePath");
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step MAC" << std::endl;
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

                auto dir = pd->getDir(access::location::cpu,
                                      access::mode::read);

                for(auto& id_i : id_index){
                    int id     = id_i.first;
                    int index  = id_i.second;

                    referenceDirections[id] = dir[index];
                }

                firstStep = false;
            }

            //Write current step to output file
            outputFile << step << " ";

            //Compute the mean angular correlation
            real mac = 0;

            auto id  = pd->getId(access::location::cpu,
                                  access::mode::read);

            auto currentDirections = this->pd->getDir(access::location::cpu,access::mode::read);

            auto groupIndex        = pg->getIndexIterator(access::location::cpu);

            fori(0,pg->getNumberParticles()){
                int id_   = id[groupIndex[i]];
                int index = id_index[id_];

                Quat qRef = referenceDirections[id_];
                Quat qCur = currentDirections[index];

                real3 EzRef = qRef.getVz();
                real3 EzCur = qCur.getVz();

                real cosTheta = dot(EzRef,EzCur);

                if(cosTheta >  1) {cosTheta =  1;}
                if(cosTheta < -1) {cosTheta = -1;}

                mac += cosTheta;
            }

            mac /= pg->getNumberParticles();

            outputFile << mac << std::endl;
        }
};

}}}}

#endif
