#pragma once

#include <string>

#include "System/ExtendedSystem.cuh"

#include "Utils/Preprocessor/DataAccess.cuh"

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

        #define VARIABLE_IMPL(NAME, name, type) \
        virtual type get##NAME(){ \
            System::log<System::CRITICAL>("[Ensemble] %s not defined for ensemble \"%s\".", \
                                          std::string(#NAME).c_str(), subType.c_str()); \
            throw std::runtime_error("Variable not defined for ensemble."); \
        } \
        virtual void set##NAME(type value){ \
            System::log<System::CRITICAL>("[Ensemble] %s not defined for ensemble \"%s\".", \
                                          std::string(#NAME).c_str(), subType.c_str()); \
        }

        // We have to use an intermediate macro to evaluate the previous expressions.
        // Otherwise, the preprocessor will not expand the __DATA_CAPS__ and __DATA_NAME__ macros.
        #define VARIABLE_EVAL(NAME, name, type) VARIABLE_IMPL(NAME, name, type)

        #define VARIABLE(r, data, seq) \
            VARIABLE_EVAL(__DATA_CAPS__(seq), __DATA_NAME__(seq), __DATA_TYPE__(seq))

        __MACRO_OVER_ENSEMBLE__(VARIABLE)

        #undef VARIABLE
        #undef VARIABLE_EVAL
        #undef VARIABLE_IMPL

        virtual void updateDataEntry(DataEntry data) = 0;


};

}}}
