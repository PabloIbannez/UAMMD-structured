#ifndef __DISTANCEBETWEENCENTERSOFMASS_MEASURE__
#define __DISTANCEBETWEENCENTERSOFMASS_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class DistanceBetweenCentersOfMass: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        std::vector<std::pair<std::shared_ptr<ParticleGroup>, std::shared_ptr<ParticleGroup>>> groups;
        std::vector<std::pair<real, real>> totalMasses;

        int nPairs;

    public:

        DistanceBetweenCentersOfMass(std::shared_ptr<ParticleGroup>     pg,
                                     std::shared_ptr<IntegratorManager> integrator,
                                     std::shared_ptr<ForceField>                ff,
                                     DataEntry& data,
                                     std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            outputFilePath = data.getParameter<std::string>("outputFilePath");

            std::vector<std::vector<int>> i_sets = data.getData<std::vector<int>>("idSet_i");
            std::vector<std::vector<int>> j_sets = data.getData<std::vector<int>>("idSet_j");

            nPairs = i_sets.size();

            for(int i = 0; i < nPairs; ++i){
                auto i_group = std::make_shared<ParticleGroup>(i_sets[i].begin(), i_sets[i].end(),
                                                               pd,"group_i_" + std::to_string(i));
                auto j_group = std::make_shared<ParticleGroup>(j_sets[i].begin(), j_sets[i].end(),
                                                               pd,"group_j_" + std::to_string(i));
                groups.push_back(std::make_pair(i_group, j_group));
            }
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step ";
                for(int i = 0; i < nPairs; ++i){
                    outputFile << "distance_" << i << " ";
                }
            }

            for(int i = 0; i < nPairs; ++i){
                totalMasses.push_back(std::make_pair(Measures::totalMass(groups[i].first),
                                                     Measures::totalMass(groups[i].second)));
            }

        }

        void applyStep(ullint step, cudaStream_t st) override{

            Box box = gd->getEnsemble()->getBox();

            outputFile << step << " ";
            for(int i = 0; i < nPairs; ++i){
                real3 com_i = Measures::centerOfMassPos(groups[i].first,totalMasses[i].first);
                real3 com_j = Measures::centerOfMassPos(groups[i].second,totalMasses[i].second);
                real3 distance = box.apply_pbc(com_i - com_j);
                outputFile << length(distance) << " ";
            }
            outputFile << std::endl;
        }
};

}}}}

#endif
