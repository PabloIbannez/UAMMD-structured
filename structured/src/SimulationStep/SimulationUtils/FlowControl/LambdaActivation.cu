#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationUtils{

class LambdaActivation: public SimulationStepBase{

      ullint lambdaValueStep;
      std::vector<real> lambdaValues;

      ullint activationLength;

  public:

      LambdaActivation(std::shared_ptr<ParticleGroup>              pg,
  	  	       std::shared_ptr<IntegratorManager>  integrator,
  	  	       std::shared_ptr<ForceField> ff,
  	  	       DataEntry& data,
  	  	       std::string name):SimulationStepBase(pg,integrator,ff,data,name){

        lambdaValueStep = data.getParameter<ullint>("lambdaValueStep");
        lambdaValues    = data.getParameter<std::vector<real>>("lambdaValues");

        System::log<System::MESSAGE>("[LambdaActivation] lambdaValueStep: %llu",lambdaValueStep);

        std::string lambdaValuesString = "";
        for(auto& lambda:lambdaValues){
          lambdaValuesString += std::to_string(lambda) + " ";
        }
        System::log<System::MESSAGE>("[LambdaActivation] lambdaValues: %s",lambdaValuesString.c_str());

      }

      void init(cudaStream_t st) override{

          //Check if intervalStep (defined in the base class) is equal to 1,
          //this simulationStep is applied every step
          if(this->intervalStep != 1){
            System::log<System::CRITICAL>("[LambdaActivation] intervalStep (%llu) must be equal to 1 for LambdaActivation.",
                this->intervalStep);
          }

          //Check if all the lambda values are between 0 and 1
          for(auto& lambda:this->lambdaValues){
            if(lambda < 0 || lambda > 1){
              System::log<System::CRITICAL>("[LambdaActivation] lambda (%f) must be between 0 and 1.",lambda);
            }
          }

          this->activationLength = this->lambdaValues.size()*this->lambdaValueStep;
      }

      void applyStep(ullint step, cudaStream_t st) override {
          ullint localStep = step - this->startStep; //This is ensure to be positive, applyStep is called only if step >= startStep

          if( localStep < this->activationLength){
            real lambda = this->lambdaValues[localStep/this->lambdaValueStep];
            if(this->gd->getEnsemble()->getLambda() != lambda){
              this->gd->getEnsemble()->setLambda(lambda);
              cudaDeviceSynchronize();
              System::log<System::MESSAGE>("[LambdaActivation] Lambda set to %f",lambda);
            }
          } else {
            if(this->gd->getEnsemble()->getLambda() != 1){
              this->gd->getEnsemble()->setLambda(1);
              cudaDeviceSynchronize();
              System::log<System::MESSAGE>("[LambdaActivation] Lambda set to 1");
            }
          }
      }

};

}}}}

REGISTER_SIMULATION_STEP(
    FlowControl,LambdaActivation,
    uammd::structured::SimulationStep::SimulationUtils::LambdaActivation
)
