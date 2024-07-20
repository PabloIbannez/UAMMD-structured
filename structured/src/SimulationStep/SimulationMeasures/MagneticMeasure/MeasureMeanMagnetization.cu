#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

#include "Utils/Measures/MeasuresBasic.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

  //Simulation step template

  class MeasureMeanMagnetization : public SimulationStepBase{

  private:

    std::string   outputFilePath;
    std::ofstream outputFile;
    int startStep;

  public:

    MeasureMeanMagnetization(std::shared_ptr<ParticleGroup>  pg,
			      std::shared_ptr<IntegratorManager> integrator,
			      std::shared_ptr<ForceField>    ff,
			      DataEntry& data,
			      std::string name):SimulationStepBase(pg,integrator,ff,data,name){

      outputFilePath = data.getParameter<std::string>("outputFilePath");
      startStep = data.getParameter<int>("startStep", 0);
    }

    void init(cudaStream_t st) override{

      bool isFileOpen = Backup::openFile(this->sys, outputFilePath, outputFile);

      //If the file did not exist, we can write the header here.
      if(!isFileOpen){
	outputFile << "Time Mx My Mz" << std::endl;
      }
    }

    void applyStep(ullint step, cudaStream_t st) override{
      if (step>=startStep){
	real3 totalMagnet = Measures::totalMagnetization(pg,st);
	real  maxMagnet   = Measures::maxMagnetization(pg,st);
	real time = step * gd->getFundamental()->getTimeStep();
	outputFile << time << " " << totalMagnet/maxMagnet << std::endl;
      }
    }
  };

}}}}

REGISTER_SIMULATION_STEP(
    MagneticMeasure,MeasureMeanMagnetization,
    uammd::structured::SimulationStep::SimulationMeasures::MeasureMeanMagnetization
)
