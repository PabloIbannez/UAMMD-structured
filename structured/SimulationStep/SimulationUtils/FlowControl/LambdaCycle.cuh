#ifndef __LAMBDA_CYCLE__
#define __LAMBDA_CYCLE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationUtils{

class LambdaCycle: public SimulationStepBase{

      ullint activationStep;
      ullint measureStep;
      ullint pauseStep;

      ullint activationLength;

      ullint cycleLength;

      enum class State{
        INACTIVE,
        ACTIVATION,
        MEASURE,
        PAUSE
      };

      State state;

      std::vector<real> lambdaValues;

      // We perform a lambda cycle every:
      // lambda=1                        ...................                               ...................
      // lambdaValues[n]              ...                                               ...
      //     .......               ...                                               ...
      //     .......            ...                                               ...
      // lambdaValues[0]     ...                                               ...
      // lambda=0                                           ...................                               ...................
      // startStep ---------><-><-><-><-><---measureStep---><----pauseStep----><-><-><-><-><---measureStep---><----pauseStep---->
      //                      |
      //                      |
      //                      L-> activationStep

  public:

      LambdaCycle(std::shared_ptr<ParticleGroup>              pg,
  	  		        std::shared_ptr<IntegratorManager>  integrator,
  	  		        std::shared_ptr<ForceField> ff,
  	  		        DataEntry& data,
  	  		        std::string name):SimulationStepBase(pg,integrator,ff,data,name){

        activationStep = data.getParameter<ullint>("activationStep");
        measureStep    = data.getParameter<ullint>("measureStep");
        pauseStep      = data.getParameter<ullint>("pauseStep");

        lambdaValues = data.getParameter<std::vector<real>>("lambdaValues");

        System::log<System::MESSAGE>("[LambdaCycle] activationStep: %llu",activationStep);
        System::log<System::MESSAGE>("[LambdaCycle] measureStep: %llu",measureStep);
        System::log<System::MESSAGE>("[LambdaCycle] pauseStep: %llu",pauseStep);

        std::string lambdaValuesString = "";
        for(auto& lambda:lambdaValues){
          lambdaValuesString += std::to_string(lambda) + " ";
        }
        System::log<System::MESSAGE>("[LambdaCycle] lambdaValues: %s",lambdaValuesString.c_str());

      }

      void init(cudaStream_t st) override{

          state = State::INACTIVE;

          //Check if intervalStep (defined in the base class) is equal to 1,
          //this simulationStep is applied every step
          if(this->intervalStep != 1){
            System::log<System::CRITICAL>("[LambdaCycle] intervalStep (%llu) must be equal to 1 for LambdaCycle.",
                this->intervalStep);
          }

          //Check if all the lambda values are between 0 and 1
          for(auto& lambda:this->lambdaValues){
            if(lambda < 0 || lambda > 1){
              System::log<System::CRITICAL>("[LambdaCycle] lambda (%f) must be between 0 and 1.",lambda);
            }
          }

          this->activationLength = this->activationStep*this->lambdaValues.size();
          this->cycleLength      = this->activationLength + this->measureStep + this->pauseStep;
      }

      void applyStep(ullint step, cudaStream_t st) override {
          ullint localStep = step - this->startStep;

          ullint cycleStep = localStep % this->cycleLength;

          if(cycleStep < this->activationLength){
            //Activation step
            ullint lambdaStep  = cycleStep % this->activationStep;
            ullint lambdaIndex = cycleStep / this->activationStep;

            real lambda = this->lambdaValues[lambdaIndex];
            if(this->gd->getEnsemble()->getLambda() != lambda or state != State::ACTIVATION){
              this->gd->getEnsemble()->setLambda(lambda);
              state = State::ACTIVATION;
              System::log<System::MESSAGE>("[LambdaCycle] Step %llu (activation): lambda = %f",step,lambda);
            }
          } else if(cycleStep < this->activationLength + this->measureStep){
            real lambda = 1.0;
            if(this->gd->getEnsemble()->getLambda() != lambda){
              this->gd->getEnsemble()->setLambda(lambda);
              state = State::MEASURE;
              System::log<System::MESSAGE>("[LambdaCycle] Step %llu (measure): lambda = %f",step,lambda);
            }
          } else {
            real lambda = 0.0;
            if(this->gd->getEnsemble()->getLambda() != lambda){
              this->gd->getEnsemble()->setLambda(lambda);
              state = State::PAUSE;
              System::log<System::MESSAGE>("[LambdaCycle] Step %llu (pause): lambda = %f",step,lambda);
            }
          }
      }

};

}}}}

#endif
