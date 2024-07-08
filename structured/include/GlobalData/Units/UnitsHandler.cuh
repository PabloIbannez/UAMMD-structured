#pragma once

#include "System/ExtendedSystem.cuh"

#include "Utils/Preprocessor/DataAccess.cuh"

#include <string>

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

        #define CONSTANT_IMPL(NAME, name, type) \
        virtual type get##NAME(){ \
            System::log<System::CRITICAL>("[Units] %s not defined for units \"%s\".", \
                                          std::string(#name).c_str(), \
                                          subType.c_str()); \
            throw std::runtime_error("Constant not defined for units."); \
        }

        #define CONSTANT_AUX(NAME, name, type) CONSTANT_IMPL(NAME, name, type)

        #define CONSTANT(r, data, tuple) \
            CONSTANT_AUX(__DATA_CAPS__(tuple), __DATA_NAME__(tuple), __DATA_TYPE__(tuple))

        __MACRO_OVER_UNITS__(CONSTANT)

        #undef CONSTANT
        #undef CONSTANT_AUX
        #undef CONSTANT_IMPL
};

}}}


