
#ifndef __UNITS_HANDLER__
#define __UNITS_HANDLER__

#include <string>

#include "InputOutput/Input/InputFormats/InputJSON.cuh"

namespace uammd{
namespace structured{
namespace Units{

class UnitsHandler{

    protected:

        std::string subType;

    public:

        UnitsHandler(DataEntry& data){
            subType = data.getSubType();
        }

        std::string getSubType() const {
            return subType;
        }

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wreturn-type"

        virtual real getBoltzmannConstant(){
            System::log<System::CRITICAL>("[Units] BoltzmannConstant not defined for units \"%s\".",
                                          subType.c_str());
        }
        virtual real getElectricConversionFactor(){
            System::log<System::CRITICAL>("[Units] ElectricConversionFactor not defined for units \"%s\".",
                                          subType.c_str());
        }



        #pragma GCC diagnostic pop



};

}}}
#endif
