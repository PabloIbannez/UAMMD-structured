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

        // We have to use an intermediate macro to evaluate the previous expressions.
        // Otherwise, the preprocessor will not expand the __DATA_CAPS__ and __DATA_NAME__ macros.
        #define CONSTANT_EVAL(NAME, name, type) CONSTANT_IMPL(NAME, name, type)

        #define CONSTANT(r, data, seq) \
            CONSTANT_EVAL(__DATA_CAPS__(seq), __DATA_NAME__(seq), __DATA_TYPE__(seq))

        __MACRO_OVER_UNITS__(CONSTANT)

        #undef CONSTANT
        #undef CONSTANT_EVAL
        #undef CONSTANT_IMPL
};

}}}


