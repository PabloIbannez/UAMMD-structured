#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/FundamentalFactory.cuh"

namespace uammd{
namespace structured{
namespace Fundamental{

    class Time: public FundamentalHandler{

        private:

            ullint currentStep;
            double timeStep;
            double simulationTime;

            int nTimeSteps=0;

        public:

            Time(DataEntry& data):FundamentalHandler(data){

                this->setCurrentStep(data.getParameter<ullint>("currentStep",0));
                this->setSimulationTime(data.getParameter<double>("simulationTime",0.0));

                if(data.isParameterAdded("timeStep")){
                    this->setTimeStep(data.getParameter<double>("timeStep"));
                } else {
                    System::log<System::WARNING>("[Time] No timeStep specified, using 0.0 as default.");
                    timeStep = 0.0;
                }
            }

            ullint getCurrentStep()    override{return currentStep;}
            double getTimeStep()       override{if(nTimeSteps>1){
                                                  System::log<System::CRITICAL>("[Time] Time step has been changed more than once. Not should be used");
                                                }
                                                return timeStep;
                                               }
            double getSimulationTime() override{return simulationTime;}

            void setCurrentStep(ullint newStep)              override{currentStep     = newStep;}
            void setTimeStep(double newTimeStep)             override{timeStep        = newTimeStep;
                                                                      nTimeSteps     += 1;}
            void setSimulationTime(double newSimulationTime) override{simulationTime = newSimulationTime;}

            void updateDataEntry(DataEntry data){
                data.setParameter("currentStep",currentStep);
                data.setParameter("timeStep",timeStep);
                data.setParameter("simulationTime",simulationTime);
            }

    };

}}}

REGISTER_FUNDAMENTAL(
    Fundamental,Time,
    uammd::structured::Fundamental::Time
)
