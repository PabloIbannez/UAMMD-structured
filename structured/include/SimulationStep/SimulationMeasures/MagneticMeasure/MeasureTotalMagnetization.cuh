//Template: Integrator

#ifndef __INTEGRATOR_MAGNETICMEASURE_MEASURETOTALMAGNETIZATION__
#define __INTEGRATOR_MAGNETICMEASURE_MEASURETOTALMAGNETIZATION__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

  //Simulation step template

  class MeasureTotalMagnetization : public SimulationStepBase{

  private:

    std::string   outputFilePath;
    std::ofstream outputFile;
    int startStep;

  public:

    MeasureTotalMagnetization(std::shared_ptr<ParticleGroup>  pg,
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
	real3 com = Measures::totalMagnetization(pg,st);
	real time = step * gd->getFundamental()->getTimeStep();
	outputFile << time << " " << com << std::endl;
      }
    }
  };

}}}}

#endif
