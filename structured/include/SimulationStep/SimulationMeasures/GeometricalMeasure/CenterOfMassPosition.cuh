//Template identified by: SimulationStep SimulationMeasures

#ifndef __CenterOfMassPosition_MEASURE__
#define __CenterOfMassPosition_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class CenterOfMassPosition: public SimulationStepBase{

        std::string   outputFilePath;
        std::ofstream outputFile;

        real totalMass;

    public:

        CenterOfMassPosition(std::shared_ptr<ParticleGroup>             pg,
                             std::shared_ptr<IntegratorManager> integrator,
                             std::shared_ptr<ForceField>                ff,
                             DataEntry& data,
                             std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            outputFilePath = data.getParameter<std::string>("outputFilePath");
        }

        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "# step CenterOfMass_X CenterOffMass_Y CenterOffMass_Z" << std::endl;
            }

            totalMass = Measures::totalMass(pg);

        }

        void applyStep(ullint step, cudaStream_t st) override{

            outputFile << step << " ";
            real3 com = Measures::centerOfMassPos(pg,totalMass);
            outputFile << com << std::endl;
        }
};

}}}}

#endif
