
#ifndef __ENSEMBLE_HANDLER__
#define __ENSEMBLE_HANDLER__

#include <string>

#include "System/ExtendedSystem.cuh"

#include "utils/Box.cuh"


namespace uammd{
namespace structured{
namespace Ensemble{

class EnsembleHandler{

    protected:

        std::string subType;

    public:

        EnsembleHandler(DataEntry& data){
            subType = data.getSubType();
        }

        std::string getSubType() const {
            return subType;
        }

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wreturn-type"

        virtual real getLambda(){
            System::log<System::CRITICAL>("[Ensemble] Lambda not defined for ensemble \"%s\".",
                                          subType.c_str());
        }
        virtual real getTemperature(){
            System::log<System::CRITICAL>("[Ensemble] Temperature not defined for ensemble \"%s\".",
                                          subType.c_str());
        }
        virtual Box getBox(){
            System::log<System::CRITICAL>("[Ensemble] Box not defined for ensemble \"%s\".",
                                          subType.c_str());
        }

        virtual void setLambda(real value){
            System::log<System::CRITICAL>("[Ensemble] Lambda not defined for ensemble \"%s\".",
                                          subType.c_str());
        }
        virtual void setTemperature(real value){
            System::log<System::CRITICAL>("[Ensemble] Temperature not defined for ensemble \"%s\".",
                                          subType.c_str());
        }
        virtual void setBox(Box value){
            System::log<System::CRITICAL>("[Ensemble] Box not defined for ensemble \"%s\".",
                                          subType.c_str());
        }


        #pragma GCC diagnostic pop

        virtual void updateDataEntry(DataEntry data) = 0;


};

}}}
#endif
