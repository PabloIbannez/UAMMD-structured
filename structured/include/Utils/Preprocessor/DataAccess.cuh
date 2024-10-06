#pragma once

#include "third_party/boost/preprocessor.hpp"
#include "third_party/boost/preprocessor/seq/for_each.hpp"

#include "Definitions/State.cuh"

#include "Definitions/Units.cuh"
#include "Definitions/Ensemble.cuh"
#include "Definitions/Types.cuh"
#include "Definitions/Fundamental.cuh"

#define __DATA_CAPS__(seq) BOOST_PP_SEQ_ELEM(0, seq)
#define __DATA_NAME__(seq) BOOST_PP_SEQ_ELEM(1, seq)
#define __DATA_TYPE__(seq) BOOST_PP_SEQ_ELEM(2, seq)

#define __MACRO_OVER_STATE__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, STATE_PROPERTIES)

#define __MACRO_OVER_UNITS__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, UNITS_PROPERTIES)
#define __MACRO_OVER_ENSEMBLE__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, ENSEMBLE_PROPERTIES)
#define __MACRO_OVER_TYPES__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, TYPES_PROPERTIES)
#define __MACRO_OVER_FUNDAMENTAL__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, FUNDAMENTAL_PROPERTIES)

