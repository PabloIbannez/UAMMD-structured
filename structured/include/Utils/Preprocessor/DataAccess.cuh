#pragma once

#include "third_party/boost/preprocessor.hpp"
#include "third_party/boost/preprocessor/seq/for_each.hpp"
#include "third_party/boost/preprocessor/tuple/elem.hpp"

#include "Definitions/State.cuh"

#include "Definitions/Units.cuh"
#include "Definitions/Ensemble.cuh"
#include "Definitions/Types.cuh"
#include "Definitions/Fundamental.cuh"

#define __DATA_CAPS__(tuple) BOOST_PP_TUPLE_ELEM(3, 0, tuple)
#define __DATA_NAME__(tuple) BOOST_PP_TUPLE_ELEM(3, 1, tuple)
#define __DATA_TYPE__(tuple) BOOST_PP_TUPLE_ELEM(3, 2, tuple)

#define __MACRO_OVER_STATE__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, STATE_PROPERTIES)

#define __MACRO_OVER_UNITS__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, UNITS_PROPERTIES)
#define __MACRO_OVER_ENSEMBLE__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, ENSEMBLE_PROPERTIES)
#define __MACRO_OVER_TYPES__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, TYPES_PROPERTIES)
#define __MACRO_OVER_FUNDAMENTAL__(macro) BOOST_PP_SEQ_FOR_EACH(macro, _, FUNDAMENTAL_PROPERTIES)

