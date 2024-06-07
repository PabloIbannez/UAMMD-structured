#ifndef __MEAN_RADIUS__
#define __MEAN_RADIUS__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class MeanRadius: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        real totalMass;

    public:

        MeanRadius(std::shared_ptr<ParticleGroup>             pg,
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
                outputFile << "# step MeanRadius" << std::endl;
            }

            totalMass = Measures::totalMass(pg);
        }

        void applyStep(ullint step, cudaStream_t st) override{

            Box box = gd->getEnsemble()->getBox();

            outputFile << step << " ";

            real3 com = Measures::centerOfMassPos(pg,totalMass);
            real  R   = Measures::meanDistance(pg,com,box);

            outputFile << R << std::endl;
        }
};

}}}}

#endif
