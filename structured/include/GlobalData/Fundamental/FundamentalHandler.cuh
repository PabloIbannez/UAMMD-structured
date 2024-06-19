
#ifndef __FUNDAMENTAL_HANDLER__
#define __FUNDAMENTAL_HANDLER__

#include <string>

#include "System/ExtendedSystem.cuh"



namespace uammd{
namespace structured{
namespace Fundamental{

class FundamentalHandler{

    protected:

        std::string subType;

    public:

        FundamentalHandler(DataEntry& data){
            subType = data.getSubType();
        }

        std::string getSubType() const {
            return subType;
        }

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wreturn-type"

        virtual double getTimeStep(){
            System::log<System::CRITICAL>("[Fundamental] TimeStep not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual ullint getCurrentStep(){
            System::log<System::CRITICAL>("[Fundamental] CurrentStep not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual double getSimulationTime(){
            System::log<System::CRITICAL>("[Fundamental] SimulationTime not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual real getEnergyThreshold(){
            System::log<System::CRITICAL>("[Fundamental] EnergyThreshold not defined for fundamental \"%s\".",
                                          subType.c_str());
        }

        virtual void setTimeStep(double value){
            System::log<System::CRITICAL>("[Fundamental] TimeStep not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual void setCurrentStep(ullint value){
            System::log<System::CRITICAL>("[Fundamental] CurrentStep not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual void setSimulationTime(double value){
            System::log<System::CRITICAL>("[Fundamental] SimulationTime not defined for fundamental \"%s\".",
                                          subType.c_str());
        }
        virtual void setEnergyThreshold(real value){
            System::log<System::CRITICAL>("[Fundamental] EnergyThreshold not defined for fundamental \"%s\".",
                                          subType.c_str());
        }


        #pragma GCC diagnostic pop

        virtual void updateDataEntry(DataEntry data) = 0;


};

}}}
#endif
